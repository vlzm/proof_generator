"""H3 induction probe, part 2 (session 10): single fixed deletion rules all
fail (h3_induction_probe.py, rules maxval/minval/pos0/posLast/head_adj) —
growing failure counts and gaps already at n = 6. Test the weaker adaptive
claim: min over ALL n choices of deleted position (not a fixed rule) of
d(delete_i(p)) satisfies d(p) <= min_i d(delete_i(p)) + (n-1) for every p.
If even this fails, no "delete one slot, recurse" induction of this exact
form can give D_n <= D_{n-1} + (n-1) with a per-permutation reduction; only
a genuinely different embedding (not "delete an array slot") could work.

Usage: python3 experiments/h3_induction_probe2.py --nmax 9
Version: h3_induction_probe2-1.0.
"""

import argparse
import itertools
import json
import os
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "exact"))
from bfs import factorials, rank_perm  # noqa: E402

VERSION = "h3_induction_probe2-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h3_induction_probe")


def load_table(n):
    path = os.path.join(ROOT, "data", "tables", f"dist_n{n}.bin")
    with open(path, "rb") as f:
        return f.read()


def reduce_delete(p, idx):
    n = len(p)
    removed_val = p[idx]
    rest = p[:idx] + p[idx + 1:]
    return tuple(v - 1 if v > removed_val else v for v in rest)


def run(n, log):
    t0 = time.time()
    fact = factorials(n)
    fact1 = factorials(n - 1)
    dist_n = load_table(n)
    dist_n1 = load_table(n - 1)
    budget = n - 1
    cnt = 0
    max_gap = -10**9
    argmax = None
    fail_count = 0
    for p in itertools.permutations(range(n)):
        cnt += 1
        dp = dist_n[rank_perm(p, fact)]
        best_dp2 = min(dist_n1[rank_perm(reduce_delete(p, i), fact1)] for i in range(n))
        gap = dp - (best_dp2 + budget)
        if gap > max_gap:
            max_gap, argmax = gap, list(p)
        if gap > 0:
            fail_count += 1
    row = {"n": n, "count": cnt, "budget_n_minus_1": budget, "seconds": round(time.time() - t0, 1),
           "max_gap": max_gap, "argmax": argmax, "fail_count": fail_count}
    status = "OK" if max_gap <= 0 else f"FAIL x{fail_count}"
    log(f"n={n}: min-over-i(delete_i) max(d - best_d' - {budget}) = {max_gap} {status} at {argmax}; "
        f"{row['seconds']} s")
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=5)
    ap.add_argument("--nmax", type=int, default=9)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    logf = open(os.path.join(OUT, "report.log"), "a")

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + "\n")
        logf.flush()

    log(f"== {VERSION} args={vars(args)}")
    rows = [run(n, log) for n in range(args.nmin, args.nmax + 1)]
    with open(os.path.join(OUT, f"report2_n{args.nmin}_{args.nmax}.json"), "w") as f:
        json.dump({"version": VERSION, "rows": rows}, f, indent=1)


if __name__ == "__main__":
    main()
