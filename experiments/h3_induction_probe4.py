"""H3 induction probe, part 4 (session 10): the "triple cut" reduction.

Sessions so far: fixed single-slot-deletion rules (probe.py) fail badly and
growing (checked against d(p) <= d(p') + (n-1) using MIN, i.e. testing one
specific candidate). Testing EXISTENCE of a working candidate among a family
(the real question for an induction proof, where the algorithm may pick
whichever reduction works) requires MAX over candidates of d(candidate),
compared against d(p) - (n-1) — a family "has a chance" for p iff some
candidate reaches that far, i.e. iff max_c d(c) + (n-1) >= d(p).

Family (per position i to delete, "triple cut"): delete array slot i
(value v = p[i]); among the remaining n-1 values (a subset of Z_n missing
v), choose a value-cut s: read them in their natural cyclic order starting
after v, but starting from offset s within that cyclic list (n-1 choices,
mirroring the double-cut of experiments/line_model.py, now for the removed
VALUE's own cyclic order rather than pi); relabel to 0..n-2 accordingly;
then choose an array rotation k (n-1 choices, head realignment) of the
resulting array. n * (n-1) * (n-1) candidates per p.

Found by hand: plain position-deletion (MAX over i only) already passes for
n <= 6, has 1 exception at n = 7, 62 at n = 8; adding rotation closes n = 7
and leaves 1 exception at n = 8; adding the value-cut closes that last one
(d(c) = 16 exactly meets the requirement for p = (3,6,0,7,4,5,2,1), n = 8).
This script re-checks the full triple-cut family exhaustively.

Usage: python3 experiments/h3_induction_probe4.py --nmax 9
Version: h3_induction_probe4-1.0.
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

VERSION = "h3_induction_probe4-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h3_induction_probe")


def load_table(n):
    path = os.path.join(ROOT, "data", "tables", f"dist_n{n}.bin")
    with open(path, "rb") as f:
        return f.read()


def best_reduction_distance(p, n, dist_n1, fact1):
    """max over the triple-cut family of d(candidate)."""
    best = 0
    for i in range(n):
        v = p[i]
        rest = p[:i] + p[i + 1:]
        base_order = [(v + 1 + t) % n for t in range(n - 1)]
        for s in range(n - 1):
            order = base_order[s:] + base_order[:s]
            label = {val: idx for idx, val in enumerate(order)}
            relabeled = tuple(label[x] for x in rest)
            for k in range(n - 1):
                rot = relabeled[k:] + relabeled[:k]
                dv = dist_n1[rank_perm(rot, fact1)]
                if dv > best:
                    best = dv
    return best


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
        best = best_reduction_distance(p, n, dist_n1, fact1)
        gap = dp - (best + budget)
        if gap > max_gap:
            max_gap, argmax = gap, list(p)
        if gap > 0:
            fail_count += 1
    row = {"n": n, "count": cnt, "budget_n_minus_1": budget, "seconds": round(time.time() - t0, 1),
           "max_gap": max_gap, "argmax": argmax, "fail_count": fail_count}
    status = "OK" if max_gap <= 0 else f"FAIL x{fail_count}"
    log(f"n={n}: triple-cut max(d - best_d' - {budget}) = {max_gap} {status} at {argmax}; "
        f"{row['seconds']} s")
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=5)
    ap.add_argument("--nmax", type=int, default=8)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    logf = open(os.path.join(OUT, "report.log"), "a")

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + "\n")
        logf.flush()

    log(f"== {VERSION} args={vars(args)}")
    rows = [run(n, log) for n in range(args.nmin, args.nmax + 1)]
    with open(os.path.join(OUT, f"report4_n{args.nmin}_{args.nmax}.json"), "w") as f:
        json.dump({"version": VERSION, "rows": rows}, f, indent=1)


if __name__ == "__main__":
    main()
