"""H13-I, session 9: independent verification of the pair-sum identity used to
justify (and then discard) the averaging approach.

Claim (reported by a prior turn in this session, not yet written up): for a
permutation pi of Z_n and an unordered pair {i, i'} with d = (i'-i) mod n and
e = (pi(i')-pi(i)) mod n (both in {1,...,n-1}, well defined up to the
d<->n-d, e<->n-e relabelling), summing inv(w) over all n^2 double cuts (q, c)
gives an exact identity:

    sum_{q,c} inv(w) = sum_{pairs {i,i'}} [ n(d+e) - 2*d*e ]

This script checks the identity exactly (integer arithmetic, brute-force
inv(w) via the same w_j formula as experiments/line_model.py and
experiments/h13_local2d.c) against all permutations for small n and random
samples for larger n.  It does NOT attempt to use the identity for a proof;
docs/notes/h13_line_model.md session 9 section explains why averaging over
all n^2 cuts cannot give the floor((n-1)^2/4) bound (already known: the
average exceeds the bound already at n=4 for some pi, e.g. identity-like
inputs), so this is bookkeeping/verification only, not a new lemma.

Usage: python3 experiments/h13_identity_check.py [--nmax 8] [--samples 200]
Version h13_identity_check-1.0.
"""
import argparse
import itertools
import random


def inv(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def total_inv_over_cuts(pi):
    n = len(pi)
    tot = 0
    for q in range(n):
        for c in range(n):
            shift = (q + 1 - c) % n
            w = [(pi[(q + 1 + j) % n] - shift) % n for j in range(n)]
            tot += inv(w)
    return tot


def identity_sum(pi):
    n = len(pi)
    tot = 0
    for i in range(n):
        for ip in range(i + 1, n):
            d = (ip - i) % n
            e = (pi[ip] - pi[i]) % n
            tot += n * (d + e) - 2 * d * e
    return tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8, help="exhaustive up to this n")
    ap.add_argument("--nmax-sampled", type=int, default=40, help="sampled up to this n")
    ap.add_argument("--samples", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    random.seed(args.seed)
    checked_exhaustive = {}
    for n in range(4, args.nmax + 1):
        cnt = 0
        for pi in itertools.permutations(range(n)):
            a = total_inv_over_cuts(list(pi))
            b = identity_sum(list(pi))
            assert a == b, (n, pi, a, b)
            cnt += 1
        checked_exhaustive[n] = cnt
        print(f"n={n}: exhaustive OK, {cnt} permutations")

    checked_sampled = {}
    for n in range(args.nmax + 1, args.nmax_sampled + 1):
        cnt = 0
        for _ in range(args.samples):
            pi = list(range(n))
            random.shuffle(pi)
            a = total_inv_over_cuts(pi)
            b = identity_sum(pi)
            assert a == b, (n, pi, a, b)
            cnt += 1
        checked_sampled[n] = cnt
        print(f"n={n}: sampled OK, {cnt} random permutations, seed={args.seed}")

    print("identity verified independently (session 9); matches the prior turn's report")


if __name__ == "__main__":
    main()
