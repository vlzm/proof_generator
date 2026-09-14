"""experiments/h13i_probe.py — session 9 probe for H13-I.

See docs/notes/h13_line_model.md S7 for context. Companion to
experiments/h13i_row_average.c (which does the exhaustive n=10, 11 check of
part 3 below in C; this script does everything else and the small-n
sanity checks the C closed form relies on).

Definitions (identical to experiments/line_model.py / docs/notes/h13_line_model.md):
  a_j = pi[(q + 1 + j) % n],  j = 0..n-1
  w_j = (a_j - s) % n
  inv(q, s) = #{(j, k) : j < k, w_j > w_k}
  I(pi) = min_{q,s} inv(q, s)          (the H13-I quantity)
  I(q)  = min_s inv(q, s)              (row minimum, q fixed)

Parts:
  1. Gradient identities (new this session):
       inv(q+1, s) - inv(q, s) = (n-1) - 2*((pi[q+1] - s) mod n)
       inv(q, s+1) - inv(q, s) = (n-1) - 2*((pi^{-1}[s] - q - 1) mod n)
     Exhaustive 4 <= n <= 7, sampled (300 random permutations) at n = 8, 9, 12.
  2. Closed form for I(q) via rho = a^{-1} and M(s) = 2*R(s) - s*(n-1)
     (see h13i_row_average.c header), checked against brute-force min_s
     inv(q, s) for every q. Exhaustive 4 <= n <= 7, sampled at n = 8, 9.
  3. The new conjecture "sum_q I(q) <= n * floor((n-1)^2/4)" checked in
     Python (exhaustive 4 <= n <= 8; n = 9..11 is done faster in the C
     program instead, see data/runs/h13i_row_average/report.md).
  4. Three explicit "constructive family" attempts to replace the true
     per-row optimum with a cheap non-adaptive rule (so that the sum bound
     would follow from ordinary pair-counting instead of an existence
     argument) -- all three are REFUTED already at small n:
       (a) diagonal family s = q + t for a single fixed t (n candidates
           (q, q+t) for a FIXED t, summed over q): no t works for every pi;
       (b) "start = 0" family (q, s) = (i-1, pi[i]) for i = 0..n-1 (so the
           row always starts with value 0): fails already on all four
           reflections of n = 4;
       (c) single-column family: min_c sum_q inv(q, c) over one fixed
           column c: fails on almost every permutation already at n = 4.
     Exhaustive 4 <= n <= 8.

Run: python3 experiments/h13i_probe.py [--nmax N]
"""
import argparse
import itertools
import math
import random


def inv_count(seq):
    n = len(seq)
    c = 0
    for j in range(n):
        for k in range(j + 1, n):
            if seq[j] > seq[k]:
                c += 1
    return c


def inv_qs(pi, q, s):
    n = len(pi)
    w = [(pi[(q + 1 + j) % n] - s) % n for j in range(n)]
    return inv_count(w)


def floor_np1_4(n):
    return ((n - 1) ** 2) // 4


def row_min_bruteforce(pi, q, n):
    return min(inv_qs(pi, q, s) for s in range(n))


def row_min_closed_form(pi, q, n):
    a = [pi[(q + 1 + j) % n] for j in range(n)]
    rho = [0] * n
    for pos, val in enumerate(a):
        rho[val] = pos
    inv_a = inv_count(a)
    R = 0
    bestM = 0
    for s in range(1, n):
        R += rho[s - 1]
        M = 2 * R - s * (n - 1)
        if M > bestM:
            bestM = M
    return inv_a - bestM


def part1_gradients(nmax):
    print("== part 1: gradient identities ==")
    random.seed(1)
    for n in range(4, nmax + 1):
        if n <= 7:
            perms = list(itertools.permutations(range(n)))
            coverage = f"exhaustive ({len(perms)} perms)"
        else:
            perms = []
            for _ in range(300):
                p = list(range(n))
                random.shuffle(p)
                perms.append(tuple(p))
            coverage = f"sampled (300 random perms)"
        bad = 0
        for pi in perms:
            pinv = [0] * n
            for i in range(n):
                pinv[pi[i]] = i
            for q in range(n):
                for s in range(n):
                    lhs_q = inv_qs(pi, (q + 1) % n, s) - inv_qs(pi, q, s)
                    rhs_q = (n - 1) - 2 * ((pi[(q + 1) % n] - s) % n)
                    if lhs_q != rhs_q:
                        bad += 1
                    lhs_s = inv_qs(pi, q, (s + 1) % n) - inv_qs(pi, q, s)
                    rhs_s = (n - 1) - 2 * ((pinv[s] - q - 1) % n)
                    if lhs_s != rhs_s:
                        bad += 1
        print(f"  n={n}: {coverage}, mismatches={bad}")


def part2_closed_form(nmax):
    print("== part 2: closed form for I(q) vs brute force ==")
    random.seed(2)
    for n in range(4, min(nmax, 9) + 1):
        if n <= 7:
            perms = list(itertools.permutations(range(n)))
            coverage = f"exhaustive ({len(perms)} perms)"
        else:
            perms = []
            for _ in range(200):
                p = list(range(n))
                random.shuffle(p)
                perms.append(tuple(p))
            coverage = "sampled (200 random perms)"
        bad = 0
        for pi in perms:
            for q in range(n):
                if row_min_closed_form(pi, q, n) != row_min_bruteforce(pi, q, n):
                    bad += 1
        print(f"  n={n}: {coverage}, mismatches={bad}")


def part3_sum_bound(nmax):
    print("== part 3: sum_q I(q) <= n * floor((n-1)^2/4) [python, small n] ==")
    for n in range(4, min(nmax, 8) + 1):
        bound = floor_np1_4(n)
        viol = 0
        worst = None
        for pi in itertools.permutations(range(n)):
            total = sum(row_min_closed_form(pi, q, n) for q in range(n))
            if total > n * bound:
                viol += 1
            if worst is None or total > worst[0]:
                worst = (total, pi)
        print(f"  n={n}: bound={bound}, n*bound={n*bound}, checked={math.factorial(n)}, "
              f"violations={viol}, worst_sum={worst[0]}, worst_pi={worst[1]}")


def part4_failed_families(nmax):
    print("== part 4: explicit non-adaptive families (all fail) ==")
    for n in range(4, min(nmax, 8) + 1):
        bound = floor_np1_4(n)

        # (a) diagonal family s = q + t, one fixed t
        universal_t = set(range(n))
        none_work = 0
        for pi in itertools.permutations(range(n)):
            oks = set()
            for t in range(n):
                total = sum(inv_qs(pi, q, (q + t) % n) for q in range(n))
                if total <= n * bound:
                    oks.add(t)
            universal_t &= oks
            if not oks:
                none_work += 1
        print(f"  n={n} (a) diagonal s=q+t: universal t for ALL pi = {universal_t}; "
              f"perms with NO working t = {none_work}/{math.factorial(n)}")

        # (b) "start = 0" family: (q,s) = (i-1, pi[i])
        worst = -1
        worst_pi = None
        viol = 0
        for pi in itertools.permutations(range(n)):
            total = 0
            for i in range(n):
                q = (i - 1) % n
                s = pi[i]
                total += inv_qs(pi, q, s)
            if total > n * bound:
                viol += 1
            if total > worst:
                worst = total
                worst_pi = pi
        print(f"  n={n} (b) start=0 family: violations={viol}/{math.factorial(n)}, "
              f"worst_sum={worst} (n*bound={n*bound}), worst_pi={worst_pi}")

        # (c) single column: min_c sum_q inv(q, c)
        worst = -1
        worst_pi = None
        viol = 0
        for pi in itertools.permutations(range(n)):
            best_col = min(sum(inv_qs(pi, q, c) for q in range(n)) for c in range(n))
            if best_col > n * bound:
                viol += 1
            if best_col > worst:
                worst = best_col
                worst_pi = pi
        print(f"  n={n} (c) single column: violations={viol}/{math.factorial(n)}, "
              f"worst_min_col_sum={worst} (n*bound={n*bound}), worst_pi={worst_pi}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=9)
    args = ap.parse_args()
    part1_gradients(max(args.nmax, 9))
    part2_closed_form(args.nmax)
    part3_sum_bound(args.nmax)
    part4_failed_families(args.nmax)
