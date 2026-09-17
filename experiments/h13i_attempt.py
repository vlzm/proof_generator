"""H13-I attempt (session 9): torical inversion minimum I(pi) <= floor((n-1)^2/4).

H13-I (docs/notes/h13_line_model.md §6): for every permutation pi of size n,
some double cut (q, c) -- rotate positions by q, rotate values by c -- gives a
line w with inv(w) <= floor((n-1)^2/4).  I(pi) = min over the n^2 cuts.

This script does NOT prove H13-I.  It records what was checked in session 9:

1. `inv_via_xor` -- an independently re-derived, exhaustively verified
   reformulation of inv(w; q, shift) as a sum over the C(n,2) point pairs of a
   single XOR of two arc-membership bits (see docs/notes/h13i_attempt.md §1).
   This is the exact formal version of the "same-name quadrants" idea named in
   h13_line_model.md §6, checked against the direct definition for all pi at
   2 <= n <= 6 and all n^2 cuts (--verify-xor).
2. `single_axis_worst` -- exhaustive worst case of I(pi) when only VALUES are
   rotated and the position order is kept fixed at the identity rotation (n
   cuts, not n^2).  Result: insufficient alone already at n = 7 (worst = 10 >
   target = 9); by the position/value symmetry of the problem, rotating only
   positions is equally insufficient.  So both degrees of freedom (q AND c)
   are necessary; no reduction of the search from n^2 to n is possible this way.
3. `greedy_two_step_worst` -- the heuristic "pick q minimising inv(x_q) first,
   then the best c for that q" -- also insufficient (n = 7: worst = 10 > 9).
   Confirms h13_line_model.md's "local optimality is not enough" (lemma C) in
   a second, independent form.
4. `structured_scan` -- I(pi) computed exactly (all n^2 cuts, O(n log n)
   inversion count per cut) on reflections (exhaustive over h), an affine
   sample, and a random sample, for n up to 50 -- SAMPLED, not VERIFIED,
   extending confidence in H13-I well past the exhaustive VERIFIED range
   (4 <= n <= 10).  Reflections match floor((n-1)^2/4) exactly at every n
   tried; affine and random inputs stay at or below it.

Usage:
  python3 experiments/h13i_attempt.py --verify-xor --nmax 6
  python3 experiments/h13i_attempt.py --single-axis --nmax 8
  python3 experiments/h13i_attempt.py --greedy --nmax 8
  python3 experiments/h13i_attempt.py --structured --ns 11,12,13,14,15,16,18,20,25,30,40,50
Output: data/runs/h13i_attempt/report.md (written by hand from the run log,
this script only prints results).  Version h13i_attempt-1.0.
"""

import argparse
import itertools
import math
import random
import time


def inv_count(seq):
    n = len(seq)
    c = 0
    for i in range(n):
        for j in range(i + 1, n):
            if seq[i] > seq[j]:
                c += 1
    return c


def inv_count_fast(seq):
    """O(n log n) inversion count via a Fenwick tree, values in [0, n)."""
    n = len(seq)
    bit = [0] * (n + 1)

    def upd(i):
        i += 1
        while i <= n:
            bit[i] += 1
            i += i & (-i)

    def qry(i):
        i += 1
        s = 0
        while i > 0:
            s += bit[i]
            i -= i & (-i)
        return s

    total = 0
    for v in reversed(seq):
        total += qry(v - 1)
        upd(v)
    return total


def target(n):
    return (n - 1) ** 2 // 4


def line_for_cut(pi, q, c):
    n = len(pi)
    shift = (q + 1 - c) % n
    positions = [(q + 1 + j) % n for j in range(n)]
    return [(pi[p] - shift) % n for p in positions]


def I_double(pi, fast=False):
    n = len(pi)
    best = None
    counter = inv_count_fast if fast else inv_count
    for q in range(n):
        positions = [(q + 1 + j) % n for j in range(n)]
        vals = [pi[p] for p in positions]
        for c in range(n):
            shift = (q + 1 - c) % n
            w = [(v - shift) % n for v in vals]
            iv = counter(w)
            if best is None or iv < best:
                best = iv
    return best


def inv_via_xor(pi, a, shift):
    """inv(w; q=a, c such that shift=(q+1-c)%n) as a sum of pairwise XORs.

    Each unordered pair of points (i1, pi[i1]), (i2, pi[i2]) defines a
    position arc of length d = (i2-i1) % n (forward from i1) and a value arc
    of length e = (v2-v1) % n (forward from v1, v1=pi[i1], v2=pi[i2]).  Let
    X = 1 iff the position-cut edge a lies on the d-arc (i.e. a+1 is inside
    (i1, i2] going forward from i1); Y = 1 iff the value shift lies on the
    e-arc similarly (edge shift-1).  The pair is inverted in w iff X xor Y.
    """
    n = len(pi)
    total = 0
    for i1 in range(n):
        v1 = pi[i1]
        for i2 in range(i1 + 1, n):
            v2 = pi[i2]
            d = (i2 - i1) % n
            e = (v2 - v1) % n
            x = 1 if (a - i1) % n < d else 0
            y = 1 if ((shift - 1) - v1) % n < e else 0
            total += x ^ y
    return total


def verify_xor(nmax):
    ok = True
    checked = 0
    for n in range(2, nmax + 1):
        for pi in itertools.permutations(range(n)):
            for q in range(n):
                for c in range(n):
                    shift = (q + 1 - c) % n
                    direct = inv_count(line_for_cut(pi, q, c))
                    via = inv_via_xor(pi, q, shift)
                    checked += 1
                    if direct != via:
                        ok = False
                        print("MISMATCH", n, pi, q, c, direct, via)
    print(f"verify_xor: nmax={nmax} checked={checked} all_match={ok}")


def single_axis_worst(nmax):
    print("n target worst_single_value_rotation_only argmax")
    for n in range(2, nmax + 1):
        worst = 0
        arg = None
        for pi in itertools.permutations(range(n)):
            best = None
            for s in range(n):
                w = [(v + s) % n for v in pi]
                iv = inv_count(w)
                if best is None or iv < best:
                    best = iv
            if best > worst:
                worst = best
                arg = pi
        print(n, target(n), worst, arg)


def greedy_two_step_worst(nmax):
    print("n target worst_greedy argmax mismatches(greedy>exact)")
    for n in range(2, nmax + 1):
        worst = 0
        arg = None
        mism = 0
        for pi in itertools.permutations(range(n)):
            best_q, best_inv_xq = None, None
            for q in range(n):
                positions = [(q + 1 + j) % n for j in range(n)]
                x = [pi[p] for p in positions]
                iv = inv_count(x)
                if best_inv_xq is None or iv < best_inv_xq:
                    best_inv_xq, best_q = iv, q
            positions = [(best_q + 1 + j) % n for j in range(n)]
            x = [pi[p] for p in positions]
            best_s = min(inv_count([(v + s) % n for v in x]) for s in range(n))
            exact = I_double(pi)
            if best_s > exact:
                mism += 1
            if best_s > worst:
                worst, arg = best_s, pi
        print(n, target(n), worst, arg, mism)


def reflection(n, h):
    return [(h - i) % n for i in range(n)]


def affine(n, a, b):
    return [(a * i + b) % n for i in range(n)]


def structured_scan(ns, affine_b_step=None, random_samples=30, seed=12345):
    rnd = random.Random(seed)
    print("n worst_reflections worst_affine worst_random target status")
    for n in ns:
        wr = max(I_double(reflection(n, h), fast=True) for h in range(n))
        wa = 0
        step = affine_b_step or max(1, n // 5)
        for a in range(1, n):
            if math.gcd(a, n) != 1:
                continue
            for b in range(0, n, step):
                v = I_double(affine(n, a, b), fast=True)
                wa = max(wa, v)
        wrand = 0
        for _ in range(random_samples):
            p = list(range(n))
            rnd.shuffle(p)
            wrand = max(wrand, I_double(p, fast=True))
        t = target(n)
        status = "OK" if max(wr, wa, wrand) <= t else "VIOLATION"
        print(n, wr, wa, wrand, t, status)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify-xor", action="store_true")
    ap.add_argument("--single-axis", action="store_true")
    ap.add_argument("--greedy", action="store_true")
    ap.add_argument("--structured", action="store_true")
    ap.add_argument("--nmax", type=int, default=6)
    ap.add_argument("--ns", type=str, default="11,12,13,14,15,16,18,20,25,30,40,50")
    args = ap.parse_args()

    t0 = time.time()
    if args.verify_xor:
        verify_xor(args.nmax)
    if args.single_axis:
        single_axis_worst(args.nmax)
    if args.greedy:
        greedy_two_step_worst(args.nmax)
    if args.structured:
        ns = [int(x) for x in args.ns.split(",") if x]
        structured_scan(ns)
    print(f"elapsed: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
