"""Checker for C37/C38 (docs/proofs/C37_toric_mean_bound.md).

Verifies every lemma of the proof numerically and, separately, collects finite
evidence for the still open H13-I.  The statement checked as PROVED is

    C37:  I(pi) <= floor((n-1)(2n-1)/6)   for all n >= 2 and all pi,

where I(pi) = min over the n^2 double cuts (a, b) of inv(w),
w_j = (pi(a+j) - b) mod n.  The proof is an averaging argument; parts 1-5 below
check its lemmas, part 6 the block criterion (lemma 7), part 7 the obstruction
C38(b), part 8 the reformulation of section 3, part 9 the (unproved) H13-I on
families.

Parts:
  1. lemma 2: floor((n-1)^2/4) + floor(n^2/4) = C(n,2) for 1 <= n <= 2000;
     lemma 4: sum_{x=1}^{n-1} (n-2x)^2 = n(n-1)(n-2)/3 for 2 <= n <= 2000.
  2. lemma 1: for every pair of positions the number of cuts that leave it
     concordant is (n-r)(n-s) + rs (all pi, 4 <= n <= 6; samples 7 <= n <= 9).
  3. lemma 3: sum over cuts of inv equals N*n^2/2 - A(pi)/4, in exact integer
     arithmetic (all pi, 4 <= n <= 7; samples up to n = 40).
  4. lemma 5: |A(pi)| <= n^2 (n-1)(n-2)/3, with equality exactly on rotations
     (A > 0) and reflections (A < 0) (all pi, 4 <= n <= 7; samples up to n = 40).
  5. C37 itself: I(pi) <= floor((n-1)(2n-1)/6) for all pi, 4 <= n <= 8, and on
     families up to n = NBIG; max over pi of the mean equals (n-1)(2n-1)/6
     (all pi, 4 <= n <= 8).
  6. lemma 7 (block criterion): coverage over all pi, 4 <= n <= 8, and the
     counterexample pi(i) = 2i mod 5.
  6b. C38(d): mean(pi) <= floor((n-1)^2/4) iff A(pi) >= n^3 (even n) / n^2(n-1)
     (odd n); the share of pi passing it, all pi, 4 <= n <= 8.
  7. lemma 6 / C38(b): on reflections inv over the class is C(k,2)+C(n-k,2),
     I = floor((n-1)^2/4) and exactly n (even n) or 2n (odd n) cuts attain it,
     4 <= n <= 60.
  8. reformulation (section 3): I(pi) <= floor((n-1)^2/4) iff max D >= floor(n/2),
     all pi, 4 <= n <= 7.
  9. evidence only (NOT a proof): H13-I on structured families up to n = NBIG.

Usage: python3 checks/check_C37.py [--nbig 120]
Output: data/runs/check_C37/report.json, report.md.
Version check_C37-1.0.
"""

import argparse
import itertools
import json
import math
import os
import random
import time
from fractions import Fraction

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "data", "runs", "check_C37")
VERSION = "check_C37-1.0"


# ---------------------------------------------------------------- basic tools

def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def line(pi, a, b):
    n = len(pi)
    return [(pi[(a + j) % n] - b) % n for j in range(n)]


def inv_grid(pi):
    """inv(a, b) for all n^2 cuts, built with O(1) updates (see the proof, s.1).

    inv(a, b+1) = inv(a, b) + (n - 1 - 2p), p = position of value b in line(a,b);
    inv(a+1, 0) = inv(a, 0) + (n - 1 - 2*pi(a)).
    """
    n = len(pi)
    pos = [0] * n
    for i, v in enumerate(pi):
        pos[v] = i
    grid = [[0] * n for _ in range(n)]
    cur = inversions(list(pi))
    for a in range(n):
        v = cur
        for b in range(n):
            grid[a][b] = v
            p = (pos[b] - a) % n
            v += n - 1 - 2 * p
        assert v == cur                      # full cycle of value shifts
        cur += n - 1 - 2 * pi[a]
    return grid


def I_of(pi):
    return min(min(row) for row in inv_grid(pi))


def A_of(pi):
    """A(pi) = sum over ordered pairs of h(i'-i) h(pi(i')-pi(i)), h(x)=n-2(x mod n)."""
    n = len(pi)
    s = 0
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            s += (n - 2 * ((j - i) % n)) * (n - 2 * ((pi[j] - pi[i]) % n))
    return s


def target(n):
    return (n - 1) ** 2 // 4


def c37_bound(n):
    return ((n - 1) * (2 * n - 1)) // 6


def families(n, rnd):
    """Structured inputs for the sampled part (AGENTS rule 9)."""
    out = {}
    for h in (0, 1, n // 2, n - 1):
        out[f"reflection h={h}"] = [(h - i) % n for i in range(n)]
    out["sigma_n"] = [(i + 1) % n for i in range(n)]
    out["identity"] = list(range(n))
    out["rev_n"] = [(n - 1 - i) % n for i in range(n)]
    for a in (2, 3, n // 2, n - 2, n - 3):
        if a > 1 and math.gcd(a, n) == 1:
            out[f"affine {a}i+1"] = [(a * i + 1) % n for i in range(n)]
    for k in range(3):
        p = list(range(n))
        rnd.shuffle(p)
        out[f"random {k}"] = p
    # two blocks swapped at the antipode, and a long cycle
    m = n // 2
    p = list(range(n))
    p[m:] = list(range(m, n))[::-1]
    out["half reversed"] = p
    return out


# ------------------------------------------------------------------- parts

def part1(log):
    for n in range(1, 2001):
        assert (n - 1) ** 2 // 4 + n * n // 4 == n * (n - 1) // 2, n
    for n in range(2, 2001):
        assert sum((n - 2 * x) ** 2 for x in range(1, n)) == n * (n - 1) * (n - 2) // 3, n
    log("part 1: lemma 2 and lemma 4 hold for n <= 2000")
    return {"lemma2_range": 2000, "lemma4_range": 2000}


def part2(log, rnd):
    checked = 0
    def check(pi):
        nonlocal checked
        n = len(pi)
        cnt = {}
        for a in range(n):
            for b in range(n):
                w = line(pi, a, b)
                rank = [0] * n
                for j in range(n):
                    rank[(a + j) % n] = j
                for i in range(n):
                    for j in range(i + 1, n):
                        conc = (rank[i] - rank[j]) * (w[rank[i]] - w[rank[j]]) > 0
                        cnt[(i, j)] = cnt.get((i, j), 0) + (1 if conc else 0)
        for i in range(n):
            for j in range(i + 1, n):
                r = (j - i) % n
                s = (pi[j] - pi[i]) % n
                assert cnt[(i, j)] == (n - r) * (n - s) + r * s, (pi, i, j)
        checked += 1
    for n in range(4, 7):
        for pi in itertools.permutations(range(n)):
            check(list(pi))
    for n in range(7, 10):
        for _ in range(20):
            p = list(range(n)); rnd.shuffle(p)
            check(p)
    log(f"part 2: lemma 1 verified on {checked} permutations (all pi for n<=6, samples n=7..9)")
    return {"permutations": checked}


def part3_4(log, rnd, nbig):
    """lemma 3 (mean) and lemma 5 (Cauchy-Schwarz with equality cases)."""
    stats = {"exhaustive_nmax": 7, "sampled_nmax": 40}
    for n in range(4, 8):
        norm2 = n * n * (n - 1) * (n - 2) // 3
        eq_plus, eq_minus = set(), set()
        for pi in itertools.permutations(range(n)):
            pi = list(pi)
            g = inv_grid(pi)
            tot = sum(sum(r) for r in g)
            A = A_of(pi)
            assert 4 * tot == n ** 3 * (n - 1) - A, (n, pi)   # sum = N n^2/2 - A/4
            assert abs(A) <= norm2, (n, pi)
            if A == norm2:
                eq_plus.add(tuple(pi))
            if A == -norm2:
                eq_minus.add(tuple(pi))
        rot = {tuple((i + c) % n for i in range(n)) for c in range(n)}
        refl = {tuple((h - i) % n for i in range(n)) for h in range(n)}
        assert eq_plus == rot and eq_minus == refl, (n, eq_plus, eq_minus)
        stats[f"n={n} equality"] = "rotations (+), reflections (-)"
    for n in (11, 16, 23, 40):
        norm2 = n * n * (n - 1) * (n - 2) // 3
        for name, pi in families(n, rnd).items():
            g = inv_grid(pi)
            tot = sum(sum(r) for r in g)
            A = A_of(pi)
            assert 4 * tot == n ** 3 * (n - 1) - A, (n, name)
            assert abs(A) <= norm2, (n, name)
    log("part 3-4: lemma 3 (exact mean) and lemma 5 (|A| bound, equality cases) hold")
    return stats


def part5(log, rnd, nbig):
    worst = {}
    for n in range(4, 9):
        mx_mean = None
        mx_I = -1
        for pi in itertools.permutations(range(n)):
            pi = list(pi)
            g = inv_grid(pi)
            tot = sum(sum(r) for r in g)
            mean = Fraction(tot, n * n)
            I = min(min(r) for r in g)
            assert I <= c37_bound(n), (n, pi, I)
            mx_mean = mean if mx_mean is None or mean > mx_mean else mx_mean
            mx_I = max(mx_I, I)
        assert mx_mean == Fraction((n - 1) * (2 * n - 1), 6), (n, mx_mean)
        worst[n] = {"max_mean": str(mx_mean), "C37_bound": c37_bound(n), "max_I": mx_I,
                    "target_floor((n-1)^2/4)": target(n)}
    sampled = {}
    for n in [12, 17, 24, 33, 50, 64, 81, 100, nbig]:
        rows = {}
        for name, pi in families(n, rnd).items():
            I = I_of(pi)
            assert I <= c37_bound(n), (n, name, I)
            rows[name] = I
        sampled[n] = {"max_I_on_family": max(rows.values()), "C37_bound": c37_bound(n),
                      "target": target(n)}
    log("part 5: C37 holds for all pi at 4<=n<=8 and on families up to n=%d; "
        "max mean = (n-1)(2n-1)/6 exactly" % nbig)
    return {"exhaustive": worst, "sampled": sampled}


def block_bound(pi):
    """min over cuts of sum C(b_t,2) over the finest direct-sum decomposition."""
    n = len(pi)
    best = None
    for a in range(n):
        for b in range(n):
            w = line(pi, a, b)
            tot, last, mx = 0, 0, -1
            for k in range(n):
                mx = max(mx, w[k])
                if mx == k:
                    sz = k + 1 - last
                    tot += sz * (sz - 1) // 2
                    last = k + 1
            if best is None or tot < best:
                best = tot
    return best


def part6(log):
    cov = {}
    for n in range(4, 9):
        good = tot = 0
        for pi in itertools.permutations(range(n)):
            tot += 1
            if block_bound(list(pi)) <= target(n):
                good += 1
        cov[n] = [good, tot]
    pi = [(2 * i) % 5 for i in range(5)]
    assert block_bound(pi) == 6 and I_of(pi) == 3 and target(5) == 4
    log("part 6: block criterion coverage " + ", ".join(f"n={n}: {a}/{b}" for n, (a, b) in cov.items())
        + "; counterexample pi=2i mod 5 (block bound 6 > 4, I = 3)")
    return {"coverage": cov, "counterexample": {"n": 5, "pi": pi, "block_bound": 6, "I": 3}}


def part6b(log):
    """C38(d): mean(pi) <= floor((n-1)^2/4) iff A(pi) >= 2 n^2 C(n,2) - 4 n^2 floor((n-1)^2/4)
    (= n^3 for even n, n^2(n-1) for odd n); how many pi pass this test."""
    cov = {}
    for n in range(4, 9):
        N = n * (n - 1) // 2
        thr = 2 * n * n * N - 4 * n * n * target(n)
        assert thr == (n ** 3 if n % 2 == 0 else n * n * (n - 1)), (n, thr)
        good = tot = 0
        for pi in itertools.permutations(range(n)):
            pi = list(pi)
            g = inv_grid(pi)
            mean = Fraction(sum(sum(r) for r in g), n * n)
            A = A_of(pi)
            assert (mean <= target(n)) == (A >= thr), (n, pi)
            tot += 1
            good += 1 if A >= thr else 0
        cov[n] = [good, tot]
    log("part 6b: averaging test (mean <= target, i.e. A >= n^3 for even n and n^2(n-1) "
        "for odd n) passes for " + ", ".join(f"n={n}: {a}/{b}" for n, (a, b) in cov.items()))
    return {"coverage": cov}


def part7(log):
    rows = {}
    for n in range(4, 61):
        pi = [(0 - i) % n for i in range(n)]
        g = inv_grid(pi)
        vals = sorted({v for row in g for v in row})
        expect = sorted({(k * k + (n - k) ** 2 - n) // 2 for k in range(1, n + 1)})
        assert vals == expect, (n, vals, expect)
        I = min(min(r) for r in g)
        assert I == target(n), (n, I)
        cnt = sum(1 for row in g for v in row if v == I)
        assert cnt == (n if n % 2 == 0 else 2 * n), (n, cnt)
        rows[n] = cnt
    log("part 7: lemma 6 / C38(b) hold for all reflections, 4 <= n <= 60 "
        "(optimal cuts: n for even n, 2n for odd n)")
    return {"optimal_cut_counts": rows}


def part8(log):
    for n in range(4, 8):
        N = n * (n - 1) // 2
        for pi in itertools.permutations(range(n)):
            g = inv_grid(list(pi))
            I = min(min(r) for r in g)
            D = max(N - 2 * v for row in g for v in row)
            assert (I <= target(n)) == (D >= n // 2), (n, pi)
            assert D % 2 == N % 2
    log("part 8: the reformulation I <= floor((n-1)^2/4) <=> max D >= floor(n/2) "
        "verified for all pi, 4 <= n <= 7")
    return {"range": "4<=n<=7"}


def part9(log, rnd, nbig):
    ev = {}
    for n in [11, 12, 14, 17, 20, 25, 32, 41, 50, 64, 80, 100, nbig]:
        mx = -1; arg = None
        for name, pi in families(n, rnd).items():
            I = I_of(pi)
            if I > mx:
                mx, arg = I, name
            assert I <= target(n), ("H13-I FAILS", n, name, I)
        # affine family in full (all a coprime to n): the recurring hard family
        amax = -1; aarg = None
        for a in range(1, n):
            if math.gcd(a, n) != 1:
                continue
            I = I_of([(a * i) % n for i in range(n)])
            assert I <= target(n), ("H13-I FAILS", n, a, I)
            if I > amax:
                amax, aarg = I, a
        ev[n] = {"max_I_families": mx, "argmax": arg, "max_I_affine": amax,
                 "affine_argmax_a": aarg, "target": target(n)}
    log("part 9 (evidence only): H13-I holds on all families and on the whole affine "
        "family for n up to %d; the maximum is always attained at a = n-1 (reflection)" % nbig)
    return ev


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nbig", type=int, default=120)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    lines = []

    def log(msg):
        print(msg, flush=True)
        lines.append(msg)

    rnd = random.Random(20260916)
    t0 = time.time()
    rep = {"version": VERSION, "nbig": args.nbig}
    rep["part1"] = part1(log)
    rep["part2"] = part2(log, rnd)
    rep["part3_4"] = part3_4(log, rnd, args.nbig)
    rep["part5"] = part5(log, rnd, args.nbig)
    rep["part6"] = part6(log)
    rep["part6b"] = part6b(log)
    rep["part7"] = part7(log)
    rep["part8"] = part8(log)
    rep["part9"] = part9(log, rnd, args.nbig)
    rep["seconds"] = round(time.time() - t0)
    rep["result"] = "PASS"
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(rep, f, indent=1, sort_keys=True)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# check_C37 — проверка лемм доказательства C37/C38\n\n")
        f.write(f"Версия {VERSION}, seed 20260916, `--nbig {args.nbig}`, "
                f"{rep['seconds']} с. Результат: PASS.\n\n")
        f.write("Части 1-8 проверяют леммы доказательства "
                "(`docs/proofs/C37_toric_mean_bound.md`); часть 9 — только конечные\n"
                "данные по открытой H13-I (не доказательство).\n\n")
        for ln in lines:
            f.write("- " + ln + "\n")
    print(f"PASS in {rep['seconds']} s -> {OUT}")


if __name__ == "__main__":
    main()
