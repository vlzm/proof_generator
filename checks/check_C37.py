"""Checker for C37 (session 9): Laplacian identity for the line-model origin
function g(a, b).

g(a, b) := inv of the line obtained from pi by cutting positions after a and
subtracting b (mod n) from every value (docs/notes/h13_line_model.md §0, §7.1);
I(pi) = min_{a,b} g(a, b) is the quantity bounded by H13-I.

Claim (Lemma D, docs/notes/h13_line_model.md §7.2): for all n >= 4, pi in S_n,
a, b in Z_n (indices a+1, b+1 taken mod n):

    g(a+1,b+1) - g(a+1,b) - g(a,b+1) + g(a,b) = 2 - 2n * [pi(a) == b]

i.e. the mixed second difference of g is +2 everywhere except at the n points
of the permutation itself, where it is 2 - 2n.  This is a structural fact
about g; it does NOT by itself prove H13-I (I(pi) <= floor((n-1)^2/4)), see
docs/notes/h13_line_model.md §7.3.

Coverage: exhaustive over all pi and all (a, b) for 4 <= n <= AMAX; random
sampling (fixed seed) for larger n as a sanity check beyond the exhaustive
range (not part of the PROVED claim, which is a theorem for all n, all pi).
Usage: python3 checks/check_C37.py [--amax 8] [--sample-n 9,10,13,20] [--trials 200]
Output: data/runs/check_C37/report.json, report.md.
"""

import argparse
import itertools
import json
import os
import random
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
from moves import CORE_VERSION  # noqa: E402

VERSION = "check_C37-1.0"
OUT = os.path.join(ROOT, "data", "runs", "check_C37")


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def grid(pi, n):
    g = [[0] * n for _ in range(n)]
    for a in range(n):
        for b in range(n):
            w = [(pi[(a + j) % n] - b) % n for j in range(n)]
            g[a][b] = inversions(w)
    return g


def laplacian_ok(pi, n):
    g = grid(pi, n)
    for a in range(n):
        for b in range(n):
            d = g[(a + 1) % n][(b + 1) % n] - g[(a + 1) % n][b] - g[a][(b + 1) % n] + g[a][b]
            expected = 2 - 2 * n if pi[a] == b else 2
            if d != expected:
                return False, (a, b, d, expected)
    return True, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amax", type=int, default=8, help="exhaustive over all pi up to this n")
    ap.add_argument("--sample-n", type=str, default="9,10,13,20")
    ap.add_argument("--trials", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    lines = [f"{VERSION} core={CORE_VERSION} amax={args.amax}"]
    rows = []
    overall_ok = True

    for n in range(4, args.amax + 1):
        t0 = time.time()
        ok_all = True
        fail = None
        count = 0
        for pi in itertools.permutations(range(n)):
            ok, info = laplacian_ok(pi, n)
            count += 1
            if not ok:
                ok_all = False
                fail = {"pi": pi, "at": info}
                break
        elapsed = time.time() - t0
        overall_ok = overall_ok and ok_all
        rows.append({"n": n, "mode": "exhaustive", "checked": count,
                      "status": "PASS" if ok_all else "FAIL", "fail": fail,
                      "elapsed_s": round(elapsed, 2)})
        lines.append(f"n={n} exhaustive ({count} pi): {'PASS' if ok_all else 'FAIL'} ({elapsed:.1f}s)")

    random.seed(args.seed)
    for n in [int(x) for x in args.sample_n.split(",") if x.strip()]:
        t0 = time.time()
        ok_all = True
        fail = None
        for _ in range(args.trials):
            pi = list(range(n))
            random.shuffle(pi)
            pi = tuple(pi)
            ok, info = laplacian_ok(pi, n)
            if not ok:
                ok_all = False
                fail = {"pi": pi, "at": info}
                break
        elapsed = time.time() - t0
        overall_ok = overall_ok and ok_all
        rows.append({"n": n, "mode": f"random({args.trials}, seed={args.seed})",
                      "checked": args.trials, "status": "PASS" if ok_all else "FAIL",
                      "fail": fail, "elapsed_s": round(elapsed, 2)})
        lines.append(f"n={n} random ({args.trials} trials): {'PASS' if ok_all else 'FAIL'} ({elapsed:.1f}s)")

    lines.append(f"OVERALL: {'PASS' if overall_ok else 'FAIL'}")
    for line in lines:
        print(line)

    report = {"version": VERSION, "core": CORE_VERSION, "claim": "C37 (Lemma D)",
               "overall": "PASS" if overall_ok else "FAIL", "rows": rows}
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(report, f, indent=2, default=list)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write(f"# check_C37 ({VERSION}, core {CORE_VERSION})\n\n")
        f.write("\n".join(lines) + "\n")

    sys.exit(0 if overall_ok else 1)


if __name__ == "__main__":
    main()
