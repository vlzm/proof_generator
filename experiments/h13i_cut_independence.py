"""h13i_cut_independence.py (h13i_cut_independence-1.0)

Session 9 attempt at H13-I (PLAN §0, docs/notes/h13_line_model.md §6):
for every pi does there exist a double cut (q, c) with inv(w) <=
floor((n-1)^2/4)?

This script checks, for small n (exhaustive, all pi), two things:

1. The "cut independence" lemma: as (q, c) range over all n^2 double cuts
   (position-cut edge q, shift c), the induced position-cut q and
   value-cut v0 = (q - c) mod n range independently over all n^2 pairs
   (q, v0) in Z_n x Z_n. Consequently, for a pair of points P=(i, pi(i)),
   P'=(j, pi(j)) with i<j, d1 = j-i, d2 = (pi(j)-pi(i)) mod n (both in
   {1,...,n-1}), the number of cuts (q, v0) for which the pair is inverted
   equals n*(d1+d2) - 2*d1*d2, independent of any other pair. This is
   verified against brute-force per-permutation counts (not just an
   average).

2. Two candidate restrictions of the n^2 cuts to only n cuts, tested for
   whether min inv over the restricted set already meets the bound:
   (a) "start value 0": for each q, take v0 = pi(q+1) (so the line starts
       with value 0);
   (b) "single rotation, value freedom only": fix q = n-1 (no position
       cut), vary v0 over all n choices.
   Both are refuted here (exhibited counterexamples): restricting to n
   cuts is not enough; the true minimum needs the joint n^2 search that
   experiments/line_profile.c already performs exhaustively for
   4 <= n <= 10 (C33/C35). This script's own reach is only 4 <= n <= 8
   (pure Python, all permutations) since it is meant to check the
   independence identity, not to extend the existing exhaustive range.

No claim here proves H13-I; see docs/notes/h13_line_model.md §7 for the
session verdict.
"""
import itertools
import argparse


def inv_of_seq(seq):
    n = len(seq)
    c = 0
    for i in range(n):
        for j in range(i + 1, n):
            if seq[i] > seq[j]:
                c += 1
    return c


def bound(n):
    return ((n - 1) ** 2) // 4


def line(pi, q, v0):
    n = len(pi)
    return [(pi[(q + 1 + j) % n] - v0) % n for j in range(n)]


def I_full(pi):
    n = len(pi)
    return min(
        inv_of_seq(line(pi, q, v0)) for q in range(n) for v0 in range(n)
    )


def pair_bad_count_formula(n, i, j, pi):
    """Number of (q, v0) in Z_n x Z_n for which pair (i, j), i<j, is inverted."""
    d1 = j - i
    d2 = (pi[j] - pi[i]) % n
    return n * (d1 + d2) - 2 * d1 * d2


def pair_bad_count_bruteforce(n, i, j, pi):
    c = 0
    for q in range(n):
        for v0 in range(n):
            w = line(pi, q, v0)
            # position of i and j in the line, and whether inverted
            pos_i = (i - q - 1) % n
            pos_j = (j - q - 1) % n
            if pos_i < pos_j:
                a, b = w[pos_i], w[pos_j]
            else:
                a, b = w[pos_j], w[pos_i]
            if a > b:
                c += 1
    return c


def check_independence(nmax):
    for n in range(4, nmax + 1):
        for pi in itertools.permutations(range(n)):
            for i in range(n):
                for j in range(i + 1, n):
                    f = pair_bad_count_formula(n, i, j, pi)
                    b = pair_bad_count_bruteforce(n, i, j, pi)
                    assert f == b, (n, pi, i, j, f, b)
        print(f"n={n}: independence formula matches brute force on all pairs, all pi")


def candidate_start_value_zero(pi):
    n = len(pi)
    best = None
    for q in range(n):
        v0 = pi[(q + 1) % n]
        iv = inv_of_seq(line(pi, q, v0))
        if best is None or iv < best:
            best = iv
    return best


def candidate_value_rotation_only(pi, qfix=0):
    n = len(pi)
    best = None
    for v0 in range(n):
        iv = inv_of_seq(line(pi, qfix, v0))
        if best is None or iv < best:
            best = iv
    return best


def scan_candidates(nmax):
    for n in range(4, nmax + 1):
        worst_a = (-1, None)
        worst_b = (-1, None)
        worst_full = (-1, None)
        for pi in itertools.permutations(range(n)):
            a = candidate_start_value_zero(pi)
            b = candidate_value_rotation_only(pi)
            f = I_full(pi)
            if a > worst_a[0]:
                worst_a = (a, pi)
            if b > worst_b[0]:
                worst_b = (b, pi)
            if f > worst_full[0]:
                worst_full = (f, pi)
        assert worst_full[0] == bound(n), (n, worst_full, bound(n))
        print(
            f"n={n} bound={bound(n)} "
            f"I_full max={worst_full[0]} (matches bound, sanity OK) "
            f"start_value_zero max={worst_a[0]} excess={worst_a[0]-bound(n)} worst={worst_a[1]} "
            f"value_rotation_only(q=0) max={worst_b[0]} excess={worst_b[0]-bound(n)} worst={worst_b[1]}"
        )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8)
    args = ap.parse_args()
    check_independence(min(args.nmax, 7))  # independence check is O(n! * n^2 * n^2); keep small
    scan_candidates(args.nmax)
