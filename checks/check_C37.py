"""Checker for C37 (session 9): the quadrant pair-inversion identity for the
double-cut ("line") model of H13-I (docs/notes/h13_line_model.md,
docs/proofs/C37_quadrant_pair_identity.md).

Claim: for n >= 2 and two points (i1, v1), (i2, v2) on the Z_n x Z_n torus
(i1 != i2, v1 != v2), with a = (i2 - i1) mod n, b = (v2 - v1) mod n, exactly
a*(n - b) + (n - a)*b of the n^2 double cuts (q, c) in Z_n x Z_n make the pair
an inversion of the line w_j = pi(q+1+j) - (q+1-c) mod n (h13_line_model.md
Section 0), where pi is any permutation with pi(i1) = v1, pi(i2) = v2 (the
count does not depend on the rest of pi, nor on the choice of pi realizing
the pair -- this is checked directly, not assumed).

Parts:
  1. Direct pair check: for 4 <= n <= NMAX, all (a, b) in {1,...,n-1}^2, at a
     handful of (i1, v1) base points (translation should not matter -- this is
     also checked, not assumed): brute-force count over all n^2 cuts against
     the formula.
  2. Whole-permutation check: for 4 <= n <= NMAX, sum over all C(n,2) pairs of
     a permutation pi of the formula equals sum_{q,c} inv(w_{q,c}) computed
     directly (cross-checks the identity's use in the note, e.g. the average
     values for id_n and sigma_n reported in h13_line_model.md Section 6).
  3. Larger n (random spot checks, no full n^2*n^2 pair enumeration): n up to
     RANDN, random (i1,v1,i2,v2) and random permutations.
  4. Lemma 2 (row-sum invariance): for 4 <= n <= PMAX, all pi, all q: the row
     sum sum_c inv(w_{q,c}) is the same for every q (proof: docs/proofs/
     C37_quadrant_pair_identity.md Lemma 2).
Usage: python3 checks/check_C37.py [--nmax 12] [--pmax 8] [--randn 60] [--trials 200]
Output: data/runs/check_C37/report.json, report.md.
"""

import argparse
import itertools
import json
import os
import random
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
VERSION = "check_C37-1.0"
OUT = os.path.join(ROOT, "data", "runs", "check_C37")


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def brute_pair_inverted_count(n, i1, v1, i2, v2):
    """Count over all n^2 cuts (q, c) whether the pair {i1, i2} (values v1,
    v2, independent of the rest of any permutation) is an inversion of the
    line."""
    cnt = 0
    for q in range(n):
        j1 = (i1 - q - 1) % n
        j2 = (i2 - q - 1) % n
        for c in range(n):
            shift = (q + 1 - c) % n
            w1 = (v1 - shift) % n
            w2 = (v2 - shift) % n
            if (j1 < j2 and w1 > w2) or (j2 < j1 and w2 > w1):
                cnt += 1
    return cnt


def formula(n, i1, v1, i2, v2):
    a = (i2 - i1) % n
    b = (v2 - v1) % n
    return a * (n - b) + (n - a) * b


def brute_sum_inv_over_all_cuts(pi):
    n = len(pi)
    total = 0
    for q in range(n):
        for c in range(n):
            shift = (q + 1 - c) % n
            w = [(pi[(q + 1 + j) % n] - shift) % n for j in range(n)]
            total += inversions(w)
    return total


def formula_sum_over_pairs(pi):
    n = len(pi)
    total = 0
    for i1 in range(n):
        for i2 in range(i1 + 1, n):
            total += formula(n, i1, pi[i1], i2, pi[i2])
    return total


def part1(nmax, log):
    ok = True
    checked = 0
    for n in range(4, nmax + 1):
        bases = [(0, 0), (1, 2), (n - 1, n - 1), (2, n - 3) if n >= 5 else (0, 1)]
        for i1, v1 in bases:
            for a in range(1, n):
                for b in range(1, n):
                    i2 = (i1 + a) % n
                    v2 = (v1 + b) % n
                    got = brute_pair_inverted_count(n, i1, v1, i2, v2)
                    want = formula(n, i1, v1, i2, v2)
                    checked += 1
                    if got != want:
                        ok = False
                        log(f"FAIL part1: n={n} i1={i1} v1={v1} a={a} b={b} "
                            f"got={got} want={want}")
    return ok, checked


def part2(pmax, log):
    ok = True
    checked = 0
    for n in range(4, pmax + 1):
        for pi in itertools.permutations(range(n)):
            got = brute_sum_inv_over_all_cuts(pi)
            want = formula_sum_over_pairs(pi)
            checked += 1
            if got != want:
                ok = False
                log(f"FAIL part2: n={n} pi={pi} got={got} want={want}")
    return ok, checked


def part4(pmax, log):
    """Lemma 2: sum_c inv(w_{q,c}) is the same for every q (fixed pi)."""
    ok = True
    checked = 0
    for n in range(4, pmax + 1):
        for pi in itertools.permutations(range(n)):
            row_sums = []
            for q in range(n):
                total = 0
                for c in range(n):
                    shift = (q + 1 - c) % n
                    w = [(pi[(q + 1 + j) % n] - shift) % n for j in range(n)]
                    total += inversions(w)
                row_sums.append(total)
            checked += 1
            if len(set(row_sums)) != 1:
                ok = False
                log(f"FAIL part4 (row-sum invariance): n={n} pi={pi} "
                    f"row_sums={row_sums}")
    return ok, checked


def part3(randn, trials, seed, log):
    rng = random.Random(seed)
    ok = True
    checked = 0
    for _ in range(trials):
        n = rng.randint(11, randn)
        i1 = rng.randrange(n)
        i2 = rng.randrange(n)
        while i2 == i1:
            i2 = rng.randrange(n)
        v1 = rng.randrange(n)
        v2 = rng.randrange(n)
        while v2 == v1:
            v2 = rng.randrange(n)
        got = brute_pair_inverted_count(n, i1, v1, i2, v2)
        want = formula(n, i1, v1, i2, v2)
        checked += 1
        if got != want:
            ok = False
            log(f"FAIL part3: n={n} i1={i1} v1={v1} i2={i2} v2={v2} "
                f"got={got} want={want}")
        # also a random full permutation, whole-sum check
        pi = list(range(n))
        rng.shuffle(pi)
        got_s = brute_sum_inv_over_all_cuts(pi)
        want_s = formula_sum_over_pairs(pi)
        checked += 1
        if got_s != want_s:
            ok = False
            log(f"FAIL part3-sum: n={n} pi={pi} got={got_s} want={want_s}")
        # row-sum invariance (Lemma 2), spot-checked at larger n too
        row_sums = []
        for q in range(n):
            total = 0
            for c in range(n):
                shift = (q + 1 - c) % n
                w = [(pi[(q + 1 + j) % n] - shift) % n for j in range(n)]
                total += inversions(w)
            row_sums.append(total)
        checked += 1
        if len(set(row_sums)) != 1:
            ok = False
            log(f"FAIL part3-rowinv: n={n} pi={pi} row_sums={row_sums}")
    return ok, checked


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=12)
    ap.add_argument("--pmax", type=int, default=8)
    ap.add_argument("--randn", type=int, default=60)
    ap.add_argument("--trials", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    lines = []

    def log(s):
        print(s)
        lines.append(s)

    t0 = time.time()
    ok1, n1 = part1(args.nmax, log)
    t1 = time.time()
    ok2, n2 = part2(args.pmax, log)
    t2 = time.time()
    ok4, n4 = part4(args.pmax, log)
    t2b = time.time()
    ok3, n3 = part3(args.randn, args.trials, args.seed, log)
    t3 = time.time()

    overall = ok1 and ok2 and ok3 and ok4
    report = {
        "version": VERSION,
        "nmax": args.nmax,
        "pmax": args.pmax,
        "randn": args.randn,
        "trials": args.trials,
        "seed": args.seed,
        "part1_pass": ok1,
        "part1_checked": n1,
        "part1_seconds": t1 - t0,
        "part2_pass": ok2,
        "part2_checked": n2,
        "part2_seconds": t2 - t1,
        "part4_pass": ok4,
        "part4_checked": n4,
        "part4_seconds": t2b - t2,
        "part3_pass": ok3,
        "part3_checked": n3,
        "part3_seconds": t3 - t2b,
        "overall_pass": overall,
        "log": lines,
    }
    with open(os.path.join(OUT, "report.json"), "w") as fh:
        json.dump(report, fh, indent=2)

    with open(os.path.join(OUT, "report.md"), "w") as fh:
        fh.write(f"# check_C37 report ({VERSION})\n\n")
        fh.write(f"nmax={args.nmax}, randn={args.randn}, trials={args.trials}, "
                 f"seed={args.seed}\n\n")
        fh.write(f"- Part 1 (pair formula vs brute force, several base points, "
                 f"all a,b, 4 <= n <= {args.nmax}): "
                 f"{'PASS' if ok1 else 'FAIL'} ({n1} checked, {t1 - t0:.2f}s)\n")
        fh.write(f"- Part 2 (whole-permutation sum vs brute force, all pi, "
                 f"4 <= n <= {args.pmax}): "
                 f"{'PASS' if ok2 else 'FAIL'} ({n2} checked, {t2 - t1:.2f}s)\n")
        fh.write(f"- Part 4 (Lemma 2, row-sum invariance, all pi, "
                 f"4 <= n <= {args.pmax}): "
                 f"{'PASS' if ok4 else 'FAIL'} ({n4} checked, {t2b - t2:.2f}s)\n")
        fh.write(f"- Part 3 (random spot checks incl. row-sum invariance, "
                 f"11 <= n <= {args.randn}): "
                 f"{'PASS' if ok3 else 'FAIL'} ({n3} checked, {t3 - t2b:.2f}s)\n\n")
        fh.write(f"Overall: {'PASS' if overall else 'FAIL'}\n")

    print(f"\nOverall: {'PASS' if overall else 'FAIL'}")
    print(f"Report: {os.path.join(OUT, 'report.md')}")
    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
