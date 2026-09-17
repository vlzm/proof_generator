"""H13-I candidate tools: exact double-cut identities and two negative results.

Two exact identities about I(pi) = min over the n^2 double cuts (a, b) of
inv(w), where w is the sequence obtained by rotating positions by a and
values by b (PLAN section 8, docs/notes/h13_line_model.md section 6):

1. Pair-cut count. For a pair of indices i != i' with position gap
   d = (i' - i) mod n and value gap e = (pi(i') - pi(i)) mod n (both in
   1..n-1), the number of cuts (a, b), out of n^2, for which the pair is
   inverted equals n*(d + e) - 2*d*e. Proof: under cut a, x'_i = (i - a) mod
   n; x'_{i'} - x'_i = d is constant, and x'_i is a bijection of a, so
   exactly n - d values of a give x'_i < x'_{i'} and d values give the
   reverse; symmetrically n - e / e for b. Since a, b vary independently,
   disagreement count = (n - d)*e + d*(n - e) = n*(d + e) - 2*d*e.

2. Telescoping in b. Fix a and let inv(a, b) be the inversion count of the
   cut-(a, b) sequence. Then inv(a, b + 1 mod n) = inv(a, b) + (n - 1) -
   2*p(a, b), where p(a, b) = (pi^{-1}(b) - a) mod n is the position (in
   the a-rotated order) of the point whose value is congruent to b. Proof:
   the cut-(a, b+1) sequence is the cut-(a, b) sequence with 1 subtracted
   mod n from every entry; this maps the unique entry equal to 0 (at
   position p) to n - 1 and decreases every other entry by exactly 1, an
   order isomorphism on the rest, so only pairs involving position p change
   order: it loses the p inversions it had as the minimum (with every
   earlier entry) and gains the n - 1 - p inversions of the new maximum
   (with every later entry).

Both are checked here against a brute-force O(n^2) cut enumeration,
exhaustively for n <= 7 and by random sampling for n = 8, 9.

Two negative results (restricted search strategies insufficient for H13-I,
i.e. worse than floor((n-1)^2/4) on reflections, with a growing gap):

A. Diagonal cuts only (a = b, n choices instead of n^2).
B. Best a, then the AVERAGE over b (rather than the best b): average_b
   inv(a, b) = inv_usual(w_a) + S(w_a) / n, where w_a(i) = pi((i + a) mod n)
   read as an actual (non-modular) sequence and S(w) = sum_i w(i) * (2*i -
   (n - 1)); minimized over a.

Both are evaluated exhaustively over all pi for n <= 7 (reflections are the
known worst case for I(pi), C33/C35).

Usage: python3 experiments/h13_cut_reformulation.py --nmax 7
       python3 experiments/h13_cut_reformulation.py --sample 8 9 --trials 300 --seed 1
Output: data/runs/h13_cut_reformulation/report.md, report.json.
Version h13_cut_reformulation-1.0.
"""

import argparse
import itertools
import json
import random
import time
from fractions import Fraction


def inv_count(w):
    n = len(w)
    c = 0
    for j in range(n):
        for k in range(j + 1, n):
            if w[j] > w[k]:
                c += 1
    return c


def inv_ab(pi, a, b):
    n = len(pi)
    xp = [(i - a) % n for i in range(n)]
    yp = [(pi[i] - b) % n for i in range(n)]
    w = [0] * n
    for i in range(n):
        w[xp[i]] = yp[i]
    return inv_count(w)


def min_inv_full(pi):
    n = len(pi)
    return min(inv_ab(pi, a, b) for a in range(n) for b in range(n))


def min_inv_diagonal(pi):
    n = len(pi)
    return min(inv_ab(pi, a, a) for a in range(n))


def avg_b_given_a(pi, a):
    n = len(pi)
    wa = [pi[(i + a) % n] for i in range(n)]
    S = sum(wa[i] * (2 * i - (n - 1)) for i in range(n))
    return inv_count(wa) + Fraction(S, n)


def min_over_a_of_avg_b(pi):
    n = len(pi)
    return min(avg_b_given_a(pi, a) for a in range(n))


def pair_count_formula(pi):
    n = len(pi)
    total = 0
    for i in range(n):
        for ip in range(i + 1, n):
            d = ip - i
            e = (pi[ip] - pi[i]) % n
            total += n * (d + e) - 2 * d * e
    return total


def pair_count_brute(pi):
    n = len(pi)
    return sum(inv_ab(pi, a, b) for a in range(n) for b in range(n))


def telescoping_ok(pi, a):
    n = len(pi)
    pi_inv = [0] * n
    for i, v in enumerate(pi):
        pi_inv[v] = i
    prev = inv_ab(pi, a, 0)
    for b in range(n):
        cur = inv_ab(pi, a, b)
        if b > 0:
            p_prev = (pi_inv[b - 1] - a) % n
            if prev + (n - 1) - 2 * p_prev != cur:
                return False
        prev = cur
    return True


def run_exhaustive(n):
    target = (n - 1) ** 2 // 4
    worst_full = worst_diag = 0
    worst_avg = Fraction(-1)
    worst_full_pi = worst_diag_pi = worst_avg_pi = None
    pair_formula_ok = True
    telescope_ok_all = True
    count = 0
    for pi in itertools.permutations(range(n)):
        pi = list(pi)
        count += 1
        mf = min_inv_full(pi)
        if mf > worst_full:
            worst_full, worst_full_pi = mf, pi
        md = min_inv_diagonal(pi)
        if md > worst_diag:
            worst_diag, worst_diag_pi = md, pi
        ma = min_over_a_of_avg_b(pi)
        if ma > worst_avg:
            worst_avg, worst_avg_pi = ma, pi
        if pair_count_formula(pi) != pair_count_brute(pi):
            pair_formula_ok = False
        for a in range(n):
            if not telescoping_ok(pi, a):
                telescope_ok_all = False
    return dict(
        n=n, target=target, coverage="exhaustive", perms_checked=count,
        worst_full=worst_full, worst_full_pi=worst_full_pi,
        worst_diag=worst_diag, worst_diag_pi=worst_diag_pi,
        worst_avg_num=worst_avg.numerator, worst_avg_den=worst_avg.denominator,
        worst_avg_pi=worst_avg_pi,
        pair_formula_ok=pair_formula_ok, telescope_ok=telescope_ok_all,
    )


def run_sampled(n, trials, seed):
    rng = random.Random(seed)
    pair_formula_ok = True
    telescope_ok_all = True
    checked = 0
    for _ in range(trials):
        pi = list(range(n))
        rng.shuffle(pi)
        checked += 1
        if pair_count_formula(pi) != pair_count_brute(pi):
            pair_formula_ok = False
        for a in range(n):
            if not telescoping_ok(pi, a):
                telescope_ok_all = False
    # also check the known worst case (reflection) explicitly
    refl = [(0 - i) % n for i in range(n)]
    target = (n - 1) ** 2 // 4
    mf = min_inv_full(refl)
    return dict(
        n=n, target=target, coverage=f"sampled ({trials} random perms, seed={seed}) + reflection",
        perms_checked=checked, reflection_min_inv_full=mf,
        pair_formula_ok=pair_formula_ok, telescope_ok=telescope_ok_all,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=7, help="exhaustive up to this n (all perms)")
    ap.add_argument("--sample", type=int, nargs="*", default=[], help="additional n values, sampled only")
    ap.add_argument("--trials", type=int, default=300)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    t0 = time.time()
    results = []
    for n in range(3, args.nmax + 1):
        results.append(run_exhaustive(n))
    for n in args.sample:
        results.append(run_sampled(n, args.trials, args.seed))
    elapsed = time.time() - t0

    for r in results:
        print(r)
    print(f"elapsed {elapsed:.1f}s")

    import os
    os.makedirs("data/runs/h13_cut_reformulation", exist_ok=True)
    with open("data/runs/h13_cut_reformulation/report.json", "w") as f:
        json.dump(dict(elapsed_s=elapsed, results=results), f, indent=2)


if __name__ == "__main__":
    main()
