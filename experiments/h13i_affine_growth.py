"""h13i_affine_growth.py -- session 9: does pair-removal induction (see
experiments/h13i_negative_results.c, test B) recover at larger n, or is its
failure at n = 8 a small-n artefact?

For every affine permutation pi(i) = a*i + b mod n (gcd(a, n) = 1), compute
  I(pi)                       -- min over q, c of inv of the relabelled line
  best_reduced = min_{p<r} I(contract(pi, p, r))
  needed       = I(pi) - best_reduced
  budget(n)    = floor((n-1)^2/4) - floor((n-3)^2/4)
and report max_{a,b} needed - budget(n) as n grows.  SAMPLED, not exhaustive
over S_n (rule 9 AGENTS.md): only the affine family, all (a, b) per n.

Usage: python3 experiments/h13i_affine_growth.py --nmax 20
Output: data/runs/h13i_negative_results/affine_growth.json.
Version h13i_affine_growth-1.0.
"""

import argparse
import json
import math
import os
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "data", "runs", "h13i_negative_results")


def inv_count(a):
    n = len(a)
    return sum(1 for i in range(n) for j in range(i + 1, n) if a[i] > a[j])


def I_of(pi):
    n = len(pi)
    best = 1 << 30
    for s in range(n):
        rot = [pi[(s + i) % n] for i in range(n)]
        for t in range(n):
            w = [(rot[i] - t) % n for i in range(n)]
            iv = inv_count(w)
            if iv < best:
                best = iv
    return best


def contract(pi, p, r):
    vals = [pi[i] for i in range(len(pi)) if i != p and i != r]
    order = {v: k for k, v in enumerate(sorted(vals))}
    return [order[v] for v in vals]


def best_pair_removal_I(pi):
    n = len(pi)
    best = 1 << 30
    for p in range(n):
        for r in range(p + 1, n):
            ip = I_of(contract(pi, p, r))
            if ip < best:
                best = ip
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=8)
    ap.add_argument("--nmax", type=int, default=20)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    rows = []
    for n in range(args.nmin, args.nmax + 1):
        m = n - 2
        budget = ((n - 1) ** 2) // 4 - ((m - 1) ** 2) // 4
        worst_needed = -1
        worst_ex = None
        checked = 0
        for a in range(1, n):
            if math.gcd(a, n) != 1:
                continue
            for b in range(n):
                pi = [(a * i + b) % n for i in range(n)]
                Ipi = I_of(pi)
                best_reduced = best_pair_removal_I(pi)
                needed = Ipi - best_reduced
                checked += 1
                if needed > worst_needed:
                    worst_needed = needed
                    worst_ex = {"a": a, "b": b, "I": Ipi, "best_reduced_I": best_reduced}
        row = {"n": n, "budget": budget, "worst_needed": worst_needed,
               "gap": worst_needed - budget, "checked_ab": checked, "argmax": worst_ex}
        rows.append(row)
        print(f"n={n} budget={budget} worst_needed={worst_needed} "
              f"gap={worst_needed - budget} argmax={worst_ex}", flush=True)
    with open(os.path.join(OUT, "affine_growth.json"), "w") as f:
        json.dump({"version": "h13i_affine_growth-1.0", "seconds": round(time.time() - t0),
                   "rows": rows}, f, indent=1)


if __name__ == "__main__":
    main()
