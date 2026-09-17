"""Checker for H13-I (session 9): I(pi) <= floor((n-1)^2/4) for all pi in S_n.

H13-I is CONJECTURED (PLAN Section 0, docs/notes/h13_line_model.md Section 6).
This checker does NOT prove it; it (a) validates the new O(n^3) exhaustive
algorithm (experiments/line_toric_max_i.c) against a fresh, independent O(n^4)
Python brute force for 4 <= n <= 8, (b) reruns/reads the exhaustive maximum of
I(pi) over all of S_n for 4 <= n <= 11 (n = 11 is new this session; n <= 10
reproduces the existing data/runs/line_profile/profile_n*.json max_I column),
and (c) runs the local-search counterexample hunt (experiments/h13i_toric_search.py)
for n = 12..14 with a modest, reproducible budget as a quick sanity re-check
(the larger one-off search from this session is logged separately, see
docs/notes/h13_line_model.md Session 9 and data/runs/h13i_verdict/).

The n = 12 EXHAUSTIVE result (max I = 30, all n = 12 reflections, ~27 min) is
a one-off run, not repeated by this checker (AGENTS rule 10); its command,
seed-free determinism and output are recorded in data/runs/h13i_verdict/report.md.

Usage: python3 checks/check_H13I.py [--n11] [--search-n 12 13 14]
Output: data/runs/check_H13I/report.json, report.md.
"""

import argparse
import itertools
import json
import os
import subprocess
import sys
import tempfile
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "experiments"))
from h13i_toric_search import run as search_run  # noqa: E402

VERSION = "check_H13I-1.0"
OUT = os.path.join(ROOT, "data", "runs", "check_H13I")


def brute_force_max_I(n):
    """Fresh, independent O(n!  * n^4) brute force (no shared code with
    experiments/line_model.py or experiments/line_toric_max_i.c): for every
    permutation, try all n^2 (q, c) cuts, recompute inv(w) from scratch."""
    best_overall = -1
    argmax = None
    for pi in itertools.permutations(range(n)):
        best = None
        for q in range(n):
            for c in range(n):
                shift = (q + 1 - c) % n
                w = [(pi[(q + 1 + j) % n] - shift) % n for j in range(n)]
                inv = sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])
                if best is None or inv < best:
                    best = inv
        if best > best_overall:
            best_overall, argmax = best, pi
    return best_overall, argmax


def run_fast_c(n, binpath, log):
    t0 = time.time()
    r = subprocess.run([binpath, str(n)], capture_output=True, text=True, timeout=1800)
    dt = time.time() - t0
    data = json.loads(r.stdout)
    log(f"n={n}: line_toric_max_i (C, O(n^3)) max_I={data['max_I']} "
        f"floor((n-1)^2/4)={data['floor_(n-1)2_4']} count={data['count']} {dt:.1f}s")
    return data, dt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n11", action="store_true", default=True, help="include the n=11 exhaustive run (~70s)")
    ap.add_argument("--no-n11", dest="n11", action="store_false")
    ap.add_argument("--search-n", type=int, nargs="*", default=[12, 13, 14])
    ap.add_argument("--search-iters", type=int, default=1500)
    ap.add_argument("--search-trials", type=int, default=4)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    lines = []

    def log(msg):
        print(msg, flush=True)
        lines.append(msg)

    log(f"== {VERSION} args={vars(args)}")

    # Compile the C program fresh, into a tmp dir (not checked into data/runs/).
    csrc = os.path.join(ROOT, "experiments", "line_toric_max_i.c")
    tmpdir = tempfile.mkdtemp(prefix="check_H13I_")
    binpath = os.path.join(tmpdir, "line_toric_max_i")
    subprocess.run(["gcc", "-O2", "-o", binpath, csrc], check=True)

    # Part A: cross-validate the O(n^3) C algorithm against an independent
    # O(n^4) Python brute force, exhaustively, for 4 <= n <= 8.
    part_a = []
    ok_a = True
    for n in range(4, 9):
        t0 = time.time()
        bf_max, bf_arg = brute_force_max_I(n)
        dt_bf = time.time() - t0
        data, dt_c = run_fast_c(n, binpath, log)
        good = (bf_max == data["max_I"] == (n - 1) ** 2 // 4)
        ok_a = ok_a and good
        part_a.append({"n": n, "brute_force_max_I": bf_max, "fast_c_max_I": data["max_I"],
                        "floor_(n-1)2_4": (n - 1) ** 2 // 4, "match": good,
                        "brute_force_seconds": round(dt_bf, 2), "fast_c_seconds": round(dt_c, 2)})
        log(f"n={n}: brute force (independent Python, O(n^4)) max_I={bf_max} "
            f"vs fast C (O(n^3)) max_I={data['max_I']} -> {'MATCH' if good else 'MISMATCH'}")

    # Part B: read existing exhaustive results for n = 9, 10 (line_profile.c,
    # session 8) and cross-check with the new O(n^3) program; add n = 11 fresh.
    part_b = []
    ok_b = True
    for n in [9, 10]:
        path = os.path.join(ROOT, "data", "runs", "line_profile", f"profile_n{n}.json")
        prev_max = None
        if os.path.exists(path):
            r = json.load(open(path))
            prev_max = max(row["I"] for row in r["by_I"])
        data, dt_c = run_fast_c(n, binpath, log)
        good = (data["max_I"] == (n - 1) ** 2 // 4) and (prev_max is None or prev_max == data["max_I"])
        ok_b = ok_b and good
        part_b.append({"n": n, "previous_max_I": prev_max, "fast_c_max_I": data["max_I"],
                        "floor_(n-1)2_4": (n - 1) ** 2 // 4, "match": good, "seconds": round(dt_c, 1)})
    if args.n11:
        data, dt_c = run_fast_c(11, binpath, log)
        good = (data["max_I"] == (11 - 1) ** 2 // 4)
        ok_b = ok_b and good
        part_b.append({"n": 11, "previous_max_I": None, "fast_c_max_I": data["max_I"],
                        "floor_(n-1)2_4": (11 - 1) ** 2 // 4, "match": good, "seconds": round(dt_c, 1),
                        "note": "new exhaustive result this session (session 9), extends C33/C36 range"})

    # Part C: local-search sanity re-check for n >= 12 (modest budget; the
    # larger one-off search is in data/runs/h13i_verdict/, not repeated here).
    part_c = []
    ok_c = True
    for n in args.search_n:
        row = search_run(n, args.search_iters, args.search_trials, seed0=2000 + n)
        ok_c = ok_c and not row["exceeds_target"]
        part_c.append(row)
        log(f"n={n}: local search (annealing, targeted at maximizing I) best={row['best']} "
            f"target={row['target']} exceeds_target={row['exceeds_target']} {row['seconds']}s")

    verdict = "PASS" if (ok_a and ok_b and ok_c) else "FAIL"
    log(f"verdict: {verdict} (this checker validates the algorithm and re-confirms no counterexample "
        f"found in the ranges tested; H13-I itself remains CONJECTURED)")

    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump({"version": VERSION, "args": vars(args), "part_a_algorithm_validation": part_a,
                   "part_b_exhaustive": part_b, "part_c_search": part_c, "verdict": verdict}, f, indent=1)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# check_H13I -- алгоритм и конечные данные (сессия 9)\n\n")
        f.write(f"Команда: `python3 checks/check_H13I.py`. Версия {VERSION}.\n\n```text\n")
        f.write("\n".join(lines) + "\n```\n")
    sys.exit(0 if verdict == "PASS" else 1)


if __name__ == "__main__":
    main()
