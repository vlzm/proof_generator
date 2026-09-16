"""h3_reduction_scan — first measurement for the H3 induction (PLAN 4.3, item 3).

The recurrence D_n <= D_{n-1} + (n-1) gives D_n <= B_n by induction, because
B_n - B_{n-1} = n - 1 and B_4 = 6.  Any induction of that shape needs a
reduction pi -> pi' from Z_n to Z_{n-1} whose simulation cost is at most n-1:

    d_n(pi) <= (n - 1) + d_{n-1}(pi')    for every pi.

This script tests the most natural reduction -- delete one element and
contract both the position circle and the value circle (the same contraction
used in C39(b) for the cut model) -- against the certified distance tables.
We use the most favourable deleted element, i.e. we test

    max_pi [ d_n(pi) - max_e d_{n-1}(pi \\ e) ] <= n - 1,

and also report the pessimistic variant with min_e.

Usage: python3 experiments/h3_reduction_scan.py [--nmin 5] [--nmax 9]
Output: data/runs/h3_reduction/report.md, report.json.
Version h3_reduction_scan-1.0.  Tables: data/tables/dist_n{n}.bin (C2v).
"""

import argparse
import itertools
import json
import os
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "exact"))
sys.path.insert(0, os.path.join(ROOT, "oracle"))

from bfs import factorials, rank_perm  # noqa: E402

VERSION = "h3_reduction_scan-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h3_reduction")


def load_table(n):
    with open(os.path.join(ROOT, "data", "tables", "dist_n%d.bin" % n), "rb") as f:
        return f.read()


def contract(perm, i):
    """Delete the element at position i, contracting positions and values."""
    v = perm[i]
    return tuple(x - 1 if x > v else x for j, x in enumerate(perm) if j != i)


def scan(n):
    big, small = load_table(n), load_table(n - 1)
    fb, fs = factorials(n), factorials(n - 1)
    best, best_arg = -10 ** 9, None
    worst, worst_arg = -10 ** 9, None
    for perm in itertools.permutations(range(n)):
        d = big[rank_perm(perm, fb)]
        red = [small[rank_perm(contract(perm, i), fs)] for i in range(n)]
        g, gm = d - max(red), d - min(red)
        if g > best:
            best, best_arg = g, perm
        if gm > worst:
            worst, worst_arg = gm, perm
    return {"n": n, "step_allowed": n - 1, "max_gap_best_element": best,
            "argmax_best": list(best_arg), "max_gap_worst_element": worst,
            "argmax_worst": list(worst_arg), "ok": best <= n - 1}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=5)
    ap.add_argument("--nmax", type=int, default=9)
    args = ap.parse_args()
    t0 = time.time()
    rows = []
    for n in range(args.nmin, args.nmax + 1):
        r = scan(n)
        rows.append(r)
        print("n=%d allowed=%d gap(best e)=%d %s gap(worst e)=%d arg=%s"
              % (n, r["step_allowed"], r["max_gap_best_element"],
                 "OK" if r["ok"] else "FAIL", r["max_gap_worst_element"],
                 r["argmax_best"]))
    rep = {"version": VERSION, "rows": rows, "seconds": round(time.time() - t0, 1),
           "result": "reduction works" if all(r["ok"] for r in rows)
                     else "reduction refuted"}
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(rep, f, indent=1, sort_keys=True)
    print(rep["result"], "%.1f s" % rep["seconds"])


if __name__ == "__main__":
    main()
