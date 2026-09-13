"""H3 feasibility probe (session 9): can a one-element reduction carry the
recurrence T_n <= T_{n-1} + (n-1) (PLAN sec. 4.3, H3)?

Reduction red_p(pi) for p in Z_n: delete position p together with its value
pi(p), re-base the remaining positions at p+1 and the remaining values at
pi(p)+1, i.e.

    red_p(pi)[i] = (pi[(p + 1 + i) mod n] - pi[p] - 1) mod n,  i = 0..n-2.

This is the only reduction that respects both cyclic orders (positions and
values) and the freedom of the target shift c.  An induction of the H3 shape
("sort the reduced permutation, then insert the deleted element") needs, for
every pi, SOME p with

    d_n(pi) <= d_{n-1}(red_p(pi)) + (n - 1),                        (*)

because B_n - B_{n-1} = n - 1 exactly (zero slack).  The probe measures
    excess(pi) = d_n(pi) - max_p d_{n-1}(red_p(pi))
against n - 1, exhaustively over S_n, from the certified tables (C2v), and
reports the histogram and the worst inputs.  It also prints the same quantity
for sigma_n and rev_n separately (the chain the induction would follow).

Usage: python3 experiments/h3_probe.py --nmax 9
Output: data/runs/h3_probe/report.json (+ hand-written report.md).
Version h3_probe-1.0; tables C2v; no reference moves are used (distances only).
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
from moves import sigma, rev  # noqa: E402

VERSION = "h3_probe-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h3_probe")


def load(n):
    with open(os.path.join(ROOT, "data", "tables", f"dist_n{n}.bin"), "rb") as fh:
        dist = fh.read()
    with open(os.path.join(ROOT, "data", "tables", f"dist_n{n}.json")) as fh:
        meta = json.load(fh)
    assert meta.get("certified"), f"table n={n} not certified"
    return dist


def reduce_at(pi, p):
    n = len(pi)
    v = pi[p]
    return [(pi[(p + 1 + i) % n] - v - 1) % n for i in range(n - 1)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=9)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    rows = []
    for n in range(5, args.nmax + 1):
        t0 = time.time()
        dn, dm = load(n), load(n - 1)
        fact, factm = factorials(n), factorials(n - 1)
        hist = {}
        worst, argw = -10 ** 9, None
        bad = 0
        for pi in itertools.permutations(range(n)):
            d = dn[rank_perm(list(pi), fact)]
            best = max(dm[rank_perm(reduce_at(pi, p), factm)] for p in range(n))
            e = d - best
            hist[e] = hist.get(e, 0) + 1
            if e > n - 1:
                bad += 1
            if e > worst:
                worst, argw = e, pi
        special = {}
        for name, pi in (("sigma", list(sigma(n))), ("rev", list(rev(n)))):
            d = dn[rank_perm(list(pi), fact)]
            vals = [dm[rank_perm(reduce_at(pi, p), factm)] for p in range(n)]
            special[name] = {"d": d, "max_red": max(vals), "min_red": min(vals),
                             "excess_max": d - max(vals), "excess_min": d - min(vals)}
        rows.append({"n": n, "need": n - 1, "max_excess": worst, "argmax": list(argw),
                     "violating": bad, "hist": dict(sorted(hist.items())),
                     "special": special, "seconds": round(time.time() - t0, 1)})
        print(f"n={n}: max excess = {worst} (нужно <= {n-1}), нарушителей {bad}, "
              f"худший {argw}; sigma {special['sigma']}, rev {special['rev']}")
    with open(os.path.join(OUT, "report.json"), "w") as fh:
        json.dump({"version": VERSION, "args": vars(args), "rows": rows}, fh, indent=1)


if __name__ == "__main__":
    main()
