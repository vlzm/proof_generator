"""H13-I contraction induction (k = 3), sampled at n = 11..20.

experiments/h13i_reduction.c checks exhaustively, for 4 <= n <= 10, whether
every pi in S_n has a k-subset of positions whose removal (contraction: drop
the k positions and their k values, close the gaps preserving order on both
circles) reduces I(pi) = min_{q,c} inv(w) by at least f(n) - f(n-k),
f(m) = floor((m-1)^2/4) -- the exact per-step increment an induction on n
would need. k = 3 passes with zero exceptions (worst case exactly meets the
budget) at every n = 4..10; k = 1 and k = 2 already fail there (see
docs/notes/h13_line_model.md sec. 7).

This script extends the k = 3 check to n = 11..20 on three families known
throughout this project to be the extremal/hard cases for I(pi) and for
related conjectures (H10, H11): all n reflections pi_h(i) = h - i mod n, all
affine maps pi(i) = a*i mod n with gcd(a, n) = 1, and 30 random permutations
per n. This is NOT exhaustive (SAMPLED, see AGENTS.md rule 9) -- it can only
refute, not confirm, the k = 3 bound in this range.

Usage: python3 experiments/h13i_reduction_sample.py --nmin 11 --nmax 20
Version h13i_reduction_sample-1.0.
"""

import argparse
import itertools
import math
import random


def I_of(pi):
    n = len(pi)
    best = None
    for q in range(n):
        order = [pi[(q + 1 + j) % n] for j in range(n)]
        for c in range(n):
            w = [(v - c) % n for v in order]
            inv = sum(1 for a in range(n) for b in range(a + 1, n) if w[a] > w[b])
            if best is None or inv < best:
                best = inv
    return best


def contract(pi, subset):
    n = len(pi)
    dropped = {pi[p] for p in subset}
    positions = [p for p in range(n) if p not in subset]
    values = [v for v in range(n) if v not in dropped]
    valmap = {v: k for k, v in enumerate(values)}
    return [valmap[pi[p]] for p in positions]


def f_bound(m):
    return ((m - 1) ** 2) // 4 if m >= 1 else 0


def best_excess_k3(pi):
    n = len(pi)
    In = I_of(pi)
    best = None
    for subset in itertools.combinations(range(n), 3):
        I2 = I_of(contract(pi, subset))
        e = In - I2
        if best is None or e < best:
            best = e
    return In, best


def reflection(n, h):
    return [(h - i) % n for i in range(n)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=11)
    ap.add_argument("--nmax", type=int, default=20)
    ap.add_argument("--nrandom", type=int, default=30)
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()

    random.seed(args.seed)
    for n in range(args.nmin, args.nmax + 1):
        budget = f_bound(n) - f_bound(n - 3)
        worst, worst_desc = None, None
        for h in range(n):
            pi = reflection(n, h)
            In, be = best_excess_k3(pi)
            if worst is None or be > worst:
                worst, worst_desc = be, f"reflection h={h}"
        for a in range(2, n):
            if math.gcd(a, n) != 1:
                continue
            pi = [(a * i) % n for i in range(n)]
            In, be = best_excess_k3(pi)
            if worst is None or be > worst:
                worst, worst_desc = be, f"affine a={a}"
        for _ in range(args.nrandom):
            pi = list(range(n))
            random.shuffle(pi)
            In, be = best_excess_k3(pi)
            if worst is None or be > worst:
                worst, worst_desc = be, f"random {pi}"
        status = "OK (within budget)" if worst <= budget else "VIOLATION"
        print(f"n={n}: budget=f(n)-f(n-3)={budget}, worst sampled excess={worst} "
              f"({worst_desc}) -> {status}")


if __name__ == "__main__":
    main()
