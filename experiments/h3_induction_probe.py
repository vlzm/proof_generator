"""H3 induction probe (session 10): test candidate reductions p (size n) ->
p' (size n-1) against the inequality d(p) <= d(p') + (n-1), using the
certified exact tables (4 <= n <= 10). This is a necessary condition for
any embedding-based induction step for D_n <= D_{n-1} + (n-1) (PLAN Sec.4.3
p.3, hypothesis H3): if a rule fails the inequality on the table, no
word-lifting construction realizing exactly that rule can work, regardless
of how clever the lifting algorithm is.

Rules (delete one array slot and close the gap; values above the removed
one are shifted down by 1 so the result is a permutation of 0..n-2):
  maxval   : delete the slot holding value n-1 (the maximum)
  minval   : delete the slot holding value 0 (the minimum)
  pos0     : delete array position 0 (whatever value is there)
  posLast  : delete array position n-1
  head_adj : delete whichever of position 0, 1 holds the larger value
             (the one "further from home" under the head-adjacent convention)

Usage: python3 experiments/h3_induction_probe.py --nmax 10
Version: h3_induction_probe-1.0.
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

VERSION = "h3_induction_probe-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h3_induction_probe")


def load_table(n):
    path = os.path.join(ROOT, "data", "tables", f"dist_n{n}.bin")
    with open(path, "rb") as f:
        return f.read()


def reduce_delete(p, idx):
    """Delete array index idx, close the gap, relabel values above the
    removed one down by 1."""
    n = len(p)
    removed_val = p[idx]
    rest = p[:idx] + p[idx + 1:]
    return tuple(v - 1 if v > removed_val else v for v in rest)


def rule_maxval(p):
    n = len(p)
    return reduce_delete(p, p.index(n - 1))


def rule_minval(p):
    return reduce_delete(p, p.index(0))


def rule_pos0(p):
    return reduce_delete(p, 0)


def rule_posLast(p):
    return reduce_delete(p, len(p) - 1)


def rule_head_adj(p):
    idx = 0 if p[0] > p[1] else 1
    return reduce_delete(p, idx)


RULES = {
    "maxval": rule_maxval,
    "minval": rule_minval,
    "pos0": rule_pos0,
    "posLast": rule_posLast,
    "head_adj": rule_head_adj,
}


def run(n, log):
    t0 = time.time()
    fact = factorials(n)
    fact1 = factorials(n - 1)
    dist_n = load_table(n)
    dist_n1 = load_table(n - 1)
    budget = n - 1
    results = {name: {"max_gap": -10**9, "argmax": None, "fail_count": 0} for name in RULES}
    cnt = 0
    for p in itertools.permutations(range(n)):
        cnt += 1
        dp = dist_n[rank_perm(p, fact)]
        for name, fn in RULES.items():
            p2 = fn(p)
            dp2 = dist_n1[rank_perm(p2, fact1)]
            gap = dp - (dp2 + budget)  # <= 0 means the rule's inequality holds
            r = results[name]
            if gap > r["max_gap"]:
                r["max_gap"] = gap
                r["argmax"] = list(p)
            if gap > 0:
                r["fail_count"] += 1
    row = {"n": n, "count": cnt, "budget_n_minus_1": budget, "seconds": round(time.time() - t0, 1),
           "rules": results}
    def fmt(name, r):
        status = "OK" if r["max_gap"] <= 0 else f"FAIL x{r['fail_count']}"
        return f"{name}: max(d-d'-{budget})={r['max_gap']} {status} at {r['argmax']}"

    summary = ", ".join(fmt(name, r) for name, r in results.items())
    log(f"n={n}: {summary}; {row['seconds']} s")
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
    with open(os.path.join(OUT, f"report_n{args.nmin}_{args.nmax}.json"), "w") as f:
        json.dump({"version": VERSION, "rows": rows}, f, indent=1)


if __name__ == "__main__":
    main()
