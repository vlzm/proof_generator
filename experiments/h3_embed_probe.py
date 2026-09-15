"""H3 probe (session 9): does naive word-replay under the obvious embedding
give the O(n-1) induction step PLAN.md section 4.3 (H3) hopes for?

Take the explicit C32 word sorting sigma_{n0} (h = 1), reverse/invert it to a
word W taking id_{n0} -> sigma_{n0}, then replay the SAME move string (same
L/R/X characters, unchanged) on id_n with n = n0 + 1 (one extra element
appended at the end of the array).  Measure d(replay result, sigma_n) via
exact BFS distance (vertex-transitivity: dist(p, q) = d(p^-1 o q)) and compare
to the H3 budget n - 1.

This is a negative/structural probe, not a construction: see
docs/notes/h13i_verdict.md section 2 for the obstruction argument (L, R
rotate the whole array, so replaying an (n0)-word on an n-array reads the
circle in a different order after one full lap; the phase drift should scale
with N_rot(W), not O(1)).

Usage: python3 experiments/h3_embed_probe.py --n0max 11
Output: data/runs/h3_embed_probe/report.json, report.md.  Version 1.0.
"""

import argparse
import json
import os
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
sys.path.insert(0, os.path.join(ROOT, "exact"))
sys.path.insert(0, os.path.join(ROOT, "constructions"))
from moves import apply_word, identity, sigma, INVERSE, CORE_VERSION  # noqa: E402
from reflection_word import reflection_word  # noqa: E402
from bfs import bfs_dict  # noqa: E402

VERSION = "h3_embed_probe-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h3_embed_probe")


def reverse_word(word):
    return "".join(INVERSE[m] for m in reversed(word))


def inverse_perm(p):
    n = len(p)
    q = [0] * n
    for i, v in enumerate(p):
        q[v] = i
    return tuple(q)


def compose(a, b):
    return tuple(a[b[i]] for i in range(len(a)))


def probe(n0, log):
    info0 = reflection_word(n0, 1)
    word_sort = info0["word"]
    assert apply_word(sigma(n0), word_sort) == identity(n0)
    word0 = reverse_word(word_sort)
    assert apply_word(identity(n0), word0) == sigma(n0)

    n = n0 + 1
    id_n = identity(n)
    res_n = apply_word(id_n, word0)
    dist_n = bfs_dict(n)
    patch = dist_n[compose(inverse_perm(res_n), sigma(n))]
    budget = n - 1
    row = {"n0": n0, "n": n, "len_word": len(word0), "N_rot_word": len(word0) - word0.count("X"),
           "patch_cost": patch, "budget_n_minus_1": budget, "within_budget": bool(patch <= budget),
           "total_len": len(word0) + patch, "B_n": n * (n - 1) // 2}
    log(f"n0={n0}: len(word)={row['len_word']} N_rot={row['N_rot_word']} patch={patch} "
        f"budget={budget} {'OK' if row['within_budget'] else 'EXCEEDS'} total={row['total_len']} B_n={row['B_n']}")
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n0max", type=int, default=9)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    logf = open(os.path.join(OUT, "report.log"), "a")

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + "\n")
        logf.flush()

    log(f"== {VERSION} {CORE_VERSION} args={vars(args)}")
    t0 = time.time()
    rows = [probe(n0, log) for n0 in range(4, args.n0max + 1)]
    exceeds = sum(1 for r in rows if not r["within_budget"])
    log(f"TOTAL: {len(rows)} probes, {exceeds} exceed the n-1 budget ({time.time() - t0:.1f} s)")
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump({"version": VERSION, "core": CORE_VERSION, "args": vars(args), "rows": rows}, f, indent=1)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write(f"# h3_embed_probe report ({VERSION}, {CORE_VERSION})\n\n")
        f.write("Naive embedding: append one element at the end of the array, replay the same\n"
                "move string that sorts sigma_n0 verbatim on the (n0+1)-array, measure the\n"
                "distance from the result to sigma_n.\n\n")
        f.write("| n0 | n | len(word) | N_rot(word) | patch cost | budget n-1 | within budget? | total | B_n |\n"
                "|---|---|---|---|---|---|---|---|---|\n")
        for r in rows:
            f.write(f"| {r['n0']} | {r['n']} | {r['len_word']} | {r['N_rot_word']} | {r['patch_cost']} | "
                    f"{r['budget_n_minus_1']} | {r['within_budget']} | {r['total_len']} | {r['B_n']} |\n")
        f.write(f"\n{exceeds}/{len(rows)} probes exceed the n-1 budget -> naive replay is not the H3 mechanism.\n")


if __name__ == "__main__":
    main()
