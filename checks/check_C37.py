#!/usr/bin/env python3
"""check_C37.py — independent checker for C37 (average bound for I(pi)).

C37 (see docs/proofs/C37_toric_average_bound.md).  For a permutation pi of Z_n
and a cut (q, c) the line is w_j = (pi(q+1+j) - (q+1-c)) mod n, j = 0..n-1;
I(pi) = min over the n^2 cuts of inv(w).  Writing a = q+1, b = q+1-c the cut is
equivalent to the pair (a, b) and w_j = (pi(a+j) - b) mod n.  Claims:

  (1) sum over the n^2 cuts of inv = n^2 * C(n,2) - S(pi), where
      S(pi) = sum over ordered pairs x != y of ((y-x) mod n) * ((pi y - pi x) mod n);
  (2) S(pi) >= n^2 (n^2 - 1) / 6, with equality exactly for the reflections
      pi(i) = h - i mod n;
  (3) I(pi) <= floor(C(n,2) - S(pi)/n^2) <= floor((2n-1)(n-1)/6).

Everything here is computed from the definitions (lines are built element by
element, inversions are counted by the double loop), independently of
experiments/toric_degree.c and experiments/line_profile.c.

Usage: python3 checks/check_C37.py [--full NMAX] [--sample K] [--seed S]
Version check_C37-1.0.
"""
from __future__ import annotations

import argparse
import itertools
import random
import sys
import time
from fractions import Fraction

VERSION = "check_C37-1.0"


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def line(pi, a, b):
    n = len(pi)
    return [(pi[(a + j) % n] - b) % n for j in range(n)]


def profile(pi):
    """inv of the line for every cut (a, b), as a flat list."""
    n = len(pi)
    return [inversions(line(pi, a, b)) for a in range(n) for b in range(n)]


def S_of(pi):
    n = len(pi)
    tot = 0
    for x in range(n):
        for y in range(n):
            if x == y:
                continue
            tot += ((y - x) % n) * ((pi[y] - pi[x]) % n)
    return tot


def is_reflection(pi):
    n = len(pi)
    h = (pi[0] + 0) % n
    return all((pi[i] + i) % n == h for i in range(n))


def binom2(n):
    return n * (n - 1) // 2


def check_one(pi, fail):
    """All three claims for a single pi.  Returns (I, S)."""
    n = len(pi)
    prof = profile(pi)
    I = min(prof)
    S = S_of(pi)
    # (1) identity
    if sum(prof) != n * n * binom2(n) - S:
        fail(f"(1) identity fails for {pi}: sum={sum(prof)}, "
             f"n^2C(n,2)-S={n * n * binom2(n) - S}")
    # (2) rearrangement bound
    if 6 * S < n * n * (n * n - 1):
        fail(f"(2) S too small for {pi}: 6S={6 * S}, n^2(n^2-1)={n * n * (n * n - 1)}")
    if (6 * S == n * n * (n * n - 1)) != is_reflection(pi):
        fail(f"(2) equality case wrong for {pi}")
    # (3) average bound, exact rational arithmetic
    avg = Fraction(sum(prof), n * n)
    if I > avg:
        fail(f"(3) I exceeds the average for {pi}: I={I}, avg={avg}")
    bound_pi = (binom2(n) * n * n - S) // (n * n)          # floor(C(n,2) - S/n^2)
    weak = (2 * n - 1) * (n - 1) // 6                      # floor((2n-1)(n-1)/6)
    if I > bound_pi:
        fail(f"(3) I > floor(avg) for {pi}: I={I}, bound={bound_pi}")
    if bound_pi > weak:
        fail(f"(3) per-pi bound exceeds the uniform bound for {pi}: "
             f"{bound_pi} > {weak}")
    # sanity: C37 is weaker than the conjecture H13-I
    if weak < (n - 1) * (n - 1) // 4:
        fail(f"C37 bound is stronger than H13-I at n={n} — arithmetic error")
    return I, S


def families(n):
    """Structured inputs: identity, rotations, reflections, affine, blocks."""
    out = [tuple(range(n)), tuple((i + 1) % n for i in range(n))]
    for h in range(n):
        out.append(tuple((h - i) % n for i in range(n)))
    for a in range(1, n):
        from math import gcd
        if gcd(a, n) == 1:
            for b in (0, 1, n // 2):
                out.append(tuple((a * i + b) % n for i in range(n)))
    # two reversed arcs (the shape of the optimal line of a reflection)
    k = n // 2
    out.append(tuple(list(range(k - 1, -1, -1)) + list(range(n - 1, k - 1, -1))))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", type=int, default=7,
                    help="exhaustive over all permutations for 4 <= n <= FULL")
    ap.add_argument("--sample", type=int, default=300,
                    help="random permutations per n in the sampled range")
    ap.add_argument("--sample-nmax", type=int, default=12)
    ap.add_argument("--large", type=str, default="20,30,40,50,60",
                    help="n for the structured spot checks")
    ap.add_argument("--seed", type=int, default=20260916)
    args = ap.parse_args()

    t0 = time.time()
    errors = []

    def fail(msg):
        errors.append(msg)
        print("FAIL:", msg)

    rnd = random.Random(args.seed)
    print(f"{VERSION}: seed {args.seed}")

    # Part 1 — exhaustive.
    for n in range(4, args.full + 1):
        maxI, minS, argminS = -1, None, []
        for pi in itertools.permutations(range(n)):
            I, S = check_one(pi, fail)
            maxI = max(maxI, I)
            if minS is None or S < minS:
                minS, argminS = S, [pi]
            elif S == minS:
                argminS.append(pi)
        refl = sorted(tuple((h - i) % n for i in range(n)) for h in range(n))
        ok_refl = sorted(argminS) == refl
        print(f"  n={n} exhaustive: max I = {maxI} (H13-I bound {(n-1)**2//4}), "
              f"C37 bound {(2*n-1)*(n-1)//6}, min S = {minS} "
              f"(= n^2(n^2-1)/6 = {n*n*(n*n-1)//6}), argmin S = the {len(argminS)} "
              f"reflections: {ok_refl}")
        if not ok_refl:
            fail(f"argmin S at n={n} is not exactly the set of reflections")
        if maxI > (n - 1) ** 2 // 4:
            fail(f"H13-I itself fails at n={n}")

    # Part 2 — random samples and families.
    for n in range(args.full + 1, args.sample_nmax + 1):
        maxI = -1
        pis = families(n) + [tuple(rnd.sample(range(n), n)) for _ in range(args.sample)]
        for pi in pis:
            I, _ = check_one(pi, fail)
            maxI = max(maxI, I)
        print(f"  n={n} sampled ({len(pis)} inputs): max I = {maxI}, "
              f"H13-I bound {(n-1)**2//4}, C37 bound {(2*n-1)*(n-1)//6}")
        if maxI > (n - 1) ** 2 // 4:
            fail(f"H13-I fails on a sampled input at n={n}")

    # Part 3 — structured spot checks at large n.
    for n in [int(x) for x in args.large.split(",") if x]:
        maxI = -1
        pis = families(n) + [tuple(rnd.sample(range(n), n)) for _ in range(3)]
        for pi in pis:
            I, _ = check_one(pi, fail)
            maxI = max(maxI, I)
        print(f"  n={n} families ({len(pis)} inputs): max I = {maxI}, "
              f"H13-I bound {(n-1)**2//4}, C37 bound {(2*n-1)*(n-1)//6}")
        if maxI > (n - 1) ** 2 // 4:
            fail(f"H13-I fails on a structured input at n={n}")

    secs = time.time() - t0
    if errors:
        print(f"FAILED with {len(errors)} errors in {secs:.1f} s")
        return 1
    print(f"PASS ({secs:.1f} s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
