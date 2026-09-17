"""check_C37.py — independent checker for the toric inversion identity (C37)
and for the finite status of H13-I / conjecture C38.

Everything here is computed from the definitions only: for each double cut
(a, s) the line w_j = (pi(a+j) - s) mod n is built explicitly and its inversions
are counted by a double loop.  No increment recursion and no result of
docs/proofs/C37_toric_identity.md is used to produce the numbers; the proved
statements are only compared against them.  Exact integer / Fraction arithmetic
throughout.

Parts:
  A  identity  inv(a,s) = c(pi) + Q(a,s)/n  with c given by the closed formula
     of Lemma 2, both of its forms                        (exhaustive 4<=n<=NA)
  B  Lemma 3:  c(pi) <= (n-1)(n-2)/6, equality iff pi is a reflection
                                            (exhaustive 4<=n<=NB, families up to NL)
  C  Theorem 5: J(r) integral and min_r J(r) <= floor((n-1)^2/4)
                                            (exhaustive 4<=n<=NB, families up to NL)
  D  Proposition 6: I(pi_h) = floor((n-1)^2/4) for every reflection   (4<=n<=ND)
  E  finite data, NOT a proof: H13-I (I(pi) <= floor((n-1)^2/4)) and conjecture
     C38 (sum_a min_s inv(a,s) <= n*floor((n-1)^2/4)) exhaustively for 4<=n<=NA,
     plus the two refuted transversal variants (anti-diagonal, shifted data cut)
     with their minimal counterexamples.

Usage: python3 checks/check_C37.py [--na 8] [--nb 8] [--nl 40] [--nd 30]
(defaults: ~40 s; --na 9 is not advisable in pure Python — use
experiments/toric_inversions.c for the exhaustive range 4 <= n <= 10)
"""
from __future__ import annotations

import argparse
import itertools
import math
import random
import sys
import time
from fractions import Fraction


# ------------------------------------------------------------------ definitions


def line(pi, a, s):
    n = len(pi)
    return [(pi[(a + j) % n] - s) % n for j in range(n)]


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def displacement_sq(w):
    return sum((j - w[j]) ** 2 for j in range(len(w)))


def cuts_table(pi):
    """[(a, s, inv, Q)] over all n^2 cuts, straight from the definitions."""
    n = len(pi)
    out = []
    for a in range(n):
        for s in range(n):
            w = line(pi, a, s)
            out.append((a, s, inversions(w), displacement_sq(w)))
    return out


def c_closed(pi):
    """c(pi) = (1/n^2) sum_{i<k} (T(k-i) - T(pi(k)-pi(i)))^2  (Lemma 2)."""
    n = len(pi)
    tot = 0
    for i in range(n):
        for k in range(i + 1, n):
            tot += ((k - i) % n - (pi[k] - pi[i]) % n) ** 2
    return Fraction(tot, n * n)


def c_data_cuts(pi):
    """c(pi) = (1/(2 n^2)) sum_k Q(k, pi(k))  (second form of Lemma 2)."""
    n = len(pi)
    tot = sum(displacement_sq(line(pi, k, pi[k])) for k in range(n))
    return Fraction(tot, 2 * n * n)


def J_values(pi):
    """J(r) = c + (1/n) sum_i D(d_i - r)^2, d_i = (i - pi(i)) mod n."""
    n = len(pi)
    c = c_closed(pi)
    out = []
    for r in range(n):
        sq = 0
        for i in range(n):
            d = (i - pi[i] - r) % n
            sq += min(d, n - d) ** 2
        out.append(c + Fraction(sq, n))
    return out


def bound(n):
    return (n - 1) ** 2 // 4


def is_reflection(pi):
    n = len(pi)
    h = pi[0] % n
    return all((pi[i] + i) % n == h for i in range(n))


def toric_I(pi):
    return min(t[2] for t in cuts_table(pi))


# ------------------------------------------------------------------ families


def family_pool(n, rng, n_random):
    pool = []
    pool += [tuple((h - i) % n for i in range(n)) for h in range(n)]          # reflections
    pool += [tuple((i + b) % n for i in range(n)) for b in range(n)]          # rotations
    for a in range(1, n):
        if math.gcd(a, n) == 1:
            for b in (0, 1, n // 2, n - 1):
                pool.append(tuple((a * i + b) % n for i in range(n)))         # affine
    for _ in range(n_random):
        p = list(range(n))
        rng.shuffle(p)
        pool.append(tuple(p))
    # reflections perturbed by one and two transpositions
    for k in (1, 2):
        for _ in range(max(4, n_random // 4)):
            h = rng.randrange(n)
            p = [(h - i) % n for i in range(n)]
            for _ in range(k):
                i, j = rng.sample(range(n), 2)
                p[i], p[j] = p[j], p[i]
            pool.append(tuple(p))
    return pool


# ------------------------------------------------------------------ parts


def part_A(na):
    print(f"[A] identity inv = c + Q/n, exhaustive 4 <= n <= {na}")
    for n in range(4, na + 1):
        t0 = time.time()
        cnt = 0
        for pi in itertools.permutations(range(n)):
            c1 = c_closed(pi)
            c2 = c_data_cuts(pi)
            assert c1 == c2, ("Lemma 2 forms differ", n, pi, c1, c2)
            for a, s, inv, Q in cuts_table(pi):
                assert inv == c1 + Fraction(Q, n), ("identity", n, pi, a, s, inv, Q, c1)
                cnt += 1
        print(f"    n={n}: {math.factorial(n)} perms, {cnt} cuts checked, "
              f"{time.time() - t0:.1f} s  OK")


def part_B(nb, nl, rng, n_random):
    print(f"[B] Lemma 3: c(pi) <= (n-1)(n-2)/6, equality iff reflection")
    for n in range(4, nb + 1):
        cmax = Fraction((n - 1) * (n - 2), 6)
        eq = 0
        for pi in itertools.permutations(range(n)):
            c = c_closed(pi)
            assert c <= cmax, ("c too large", n, pi, c, cmax)
            if c == cmax:
                assert is_reflection(pi), ("equality off a reflection", n, pi)
                eq += 1
        assert eq == n, (n, eq)
        print(f"    n={n}: exhaustive, equality on exactly {eq} permutations "
              f"(= the n reflections)  OK")
    for n in range(nb + 1, nl + 1):
        cmax = Fraction((n - 1) * (n - 2), 6)
        pool = family_pool(n, rng, n_random)
        eq = 0
        for pi in pool:
            c = c_closed(pi)
            assert c <= cmax, ("c too large", n, pi, c, cmax)
            if c == cmax:
                assert is_reflection(pi), ("equality off a reflection", n, pi)
                eq += 1
        assert eq >= n
    print(f"    n={nb+1}..{nl}: families + random samples  OK")


def part_C(nb, nl, rng, n_random):
    print(f"[C] Theorem 5: J(r) integral and min_r J(r) <= floor((n-1)^2/4)")
    for n in range(4, nb + 1):
        b = bound(n)
        for pi in itertools.permutations(range(n)):
            Js = J_values(pi)
            for J in Js:
                assert J.denominator == 1, ("J not integral", n, pi, J)
            assert min(Js) <= b, ("Theorem 5 violated", n, pi, min(Js), b)
            # Lemma 4(a): J(r) is a lower bound for every cut on the diagonal r
            if n <= 7:
                for a, s, inv, Q in cuts_table(pi):
                    assert inv >= Js[(a - s) % n], ("J not a lower bound", n, pi, a, s)
        print(f"    n={n}: exhaustive"
              + (" (with the per-cut lower bound)" if n <= 7 else "") + "  OK")
    for n in range(nb + 1, nl + 1):
        b = bound(n)
        for pi in family_pool(n, rng, n_random):
            Js = J_values(pi)
            assert all(J.denominator == 1 for J in Js), (n, pi)
            assert min(Js) <= b, ("Theorem 5 violated", n, pi, min(Js), b)
    print(f"    n={nb+1}..{nl}: families + random samples  OK")


def part_D(nd):
    print(f"[D] Proposition 6: I(pi_h) = floor((n-1)^2/4) for all reflections")
    for n in range(4, nd + 1):
        b = bound(n)
        for h in range(n):
            pi = tuple((h - i) % n for i in range(n))
            I = toric_I(pi)
            assert I == b, ("reflection value", n, h, I, b)
    print(f"    4 <= n <= {nd}, all h  OK")


def part_E(na):
    print(f"[E] finite data (NOT a proof): H13-I and conjecture C38, "
          f"exhaustive 4 <= n <= {na}")
    for n in range(4, na + 1):
        b = bound(n)
        maxI, argmaxI, n_argmax, n_refl = -1, None, 0, 0
        maxrow = -1
        anti_bad, shift_bad = None, None
        for pi in itertools.permutations(range(n)):
            tab = {}
            for a, s, inv, Q in cuts_table(pi):
                tab[(a, s)] = inv
            I = min(tab.values())
            assert I <= b, ("H13-I violated", n, pi, I, b)
            if I > maxI:
                maxI, argmaxI, n_argmax, n_refl = I, pi, 0, 0
            if I == maxI:
                n_argmax += 1
                n_refl += 1 if is_reflection(pi) else 0
            rows = sum(min(tab[(a, s)] for s in range(n)) for a in range(n))
            assert rows <= n * b, ("C38 violated", n, pi, rows, n * b)
            maxrow = max(maxrow, rows)
            if anti_bad is None:
                v = min(sum(tab[(a, (t - a) % n)] for a in range(n)) for t in range(n))
                if v > n * b:
                    anti_bad = (pi, v, n * b)
            if shift_bad is None:
                v = min(sum(tab[(a, (pi[a] + g) % n)] for a in range(n)) for g in range(n))
                if v > n * b:
                    shift_bad = (pi, v, n * b)
        assert maxI == b and n_argmax == n and n_refl == n, (n, maxI, n_argmax, n_refl)
        assert maxrow == n * b, (n, maxrow)
        print(f"    n={n}: max I = {maxI} = floor((n-1)^2/4), argmax = the {n} "
              f"reflections; max sum_a min_s inv = {maxrow} = n*bound")
        print(f"        anti-diagonal average: first counterexample "
              f"{anti_bad if anti_bad else 'none'}")
        print(f"        shifted data-cut average: first counterexample "
              f"{shift_bad if shift_bad else 'none'}")


def part_F():
    print("[F] recorded counterexamples of the verdict note")
    # Q0 (min Q <= sum_v D(v)^2 - n/2 for even n) fails at n=8, pi(i)=5i
    n = 8
    pi = tuple((5 * i) % n for i in range(n))
    minQ = min(t[3] for t in cuts_table(pi))
    thr = sum(min(v, n - v) ** 2 for v in range(n)) - n // 2
    assert minQ == 56 and thr == 40, (minQ, thr)
    print(f"    n=8, pi(i)=5i mod 8: min Q = {minQ} > {thr} = sum_v D(v)^2 - n/2  "
          f"(relaxation 'drop c' is false)")
    # per-anchor form of C38 fails already at n=4
    n, tau = 4, (0, 2, 3, 1)
    V = sum((j - tau[j]) ** 2 for j in range(n))
    W = min(sum((j - ((tau[j] - sg) % n)) ** 2 for j in range(n)) for sg in range(n))
    assert 2 * W + V > 2 * n * bound(n), (W, V)
    print(f"    n=4, tau={tau}: min_sigma W + V/2 = {W + Fraction(V, 2)} > "
          f"{n * bound(n)} = n*bound (per-anchor form of C38 is false)")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--na", type=int, default=8, help="exhaustive range for parts A, E")
    ap.add_argument("--nb", type=int, default=8, help="exhaustive range for parts B, C")
    ap.add_argument("--nl", type=int, default=40, help="family range for parts B, C")
    ap.add_argument("--nd", type=int, default=30, help="range for part D")
    ap.add_argument("--random", type=int, default=20, help="random perms per n")
    ap.add_argument("--seed", type=int, default=20260917)
    args = ap.parse_args(argv)
    rng = random.Random(args.seed)
    t0 = time.time()
    part_A(args.na)
    part_B(args.nb, args.nl, rng, args.random)
    part_C(args.nb, args.nl, rng, args.random)
    part_D(args.nd)
    part_E(args.na)
    part_F()
    print(f"PASS  ({time.time() - t0:.1f} s, seed {args.seed})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
