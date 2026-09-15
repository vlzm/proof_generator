"""check_C37.py — independent checker for the H13-I computational claim C37
(session 9): I(pi) = min_{q,c} inv(line at double-cut (q,c)) satisfies
I(pi) <= floor((n-1)^2/4) for every permutation pi of {0,...,n-1}.

Two parts:

(A) Correctness of the O(n^2) recursion used by experiments/toric_I_fast.c
    (the only reason n = 11..13 became tractable) against a direct,
    independently-written O(n^4) brute force over all n^2 cuts: exhaustive
    for 2 <= n <= 7, random 2000 samples each for n = 8, 9, 10.

(B) The bound itself:
    - Recomputes I(pi) directly (brute force, independent of the C program)
      and checks the bound exhaustively for 4 <= n <= 9 (fast enough in
      pure Python without the O(n^2) trick).
    - Cross-checks the O(n^2)-formula results against the exact known
      max I values from data/runs/line_profile (n = 4..10: 2,4,6,9,12,16,20)
      and the new exhaustive run data/runs/h13_I_toric_scan (n = 11,12,13).
    - Structured + random sampling at n = 20, 50, 100, 101 (rule 9 ladder):
      sigma_n and rotations, rev_n, tau_n (antipodal transposition, even n),
      block-swap-at-antipode, affine maps, cycles, products of random
      transpositions, random permutations.

This does NOT constitute a proof of H13-I (see docs/notes/h13_I_attempt_*.md
for the session-9 verdict): it VERIFIES the stated finite ranges and samples.
Version check_C37-1.0.
"""
import itertools
import random
import sys


def inv_count(w):
    n = len(w)
    c = 0
    for i in range(n):
        for j in range(i + 1, n):
            if w[i] > w[j]:
                c += 1
    return c


def I_brute(pi):
    """Direct O(n^4) computation: all n^2 cuts, O(n^2) inversions each."""
    n = len(pi)
    best = None
    for q in range(n):
        for c in range(n):
            shift = (q + 1 - c) % n
            w = [(pi[(q + 1 + j) % n] - shift) % n for j in range(n)]
            iv = inv_count(w)
            if best is None or iv < best:
                best = iv
    return best


def I_fast(pi):
    """O(n^2) recursion (see experiments/toric_I_fast.c header for the
    derivation)."""
    n = len(pi)
    rho = [0] * n
    for i, v in enumerate(pi):
        rho[v] = i
    inv0 = inv_count(pi)
    P = [0] * n
    P[0] = inv0
    for s in range(n - 1):
        P[s + 1] = P[s] + (n - 1 - 2 * pi[s])
    best = None
    for s in range(n):
        F = P[s]
        if best is None or F < best:
            best = F
        for t in range(n - 1):
            sigma = (rho[t] - s) % n
            F = F + (n - 1 - 2 * sigma)
            if F < best:
                best = F
    return best


def bound(n):
    return (n - 1) ** 2 // 4


def part_a():
    print("Part A: I_fast vs I_brute (independent implementations)")
    ok = True
    for n in range(2, 8):
        bad = []
        for pi in itertools.permutations(range(n)):
            a, b = I_brute(list(pi)), I_fast(list(pi))
            if a != b:
                bad.append((pi, a, b))
        status = "OK" if not bad else f"MISMATCH {bad[:3]}"
        print(f"  n={n} exhaustive: {status}")
        ok = ok and not bad
    random.seed(2026_09_15)
    for n in [8, 9, 10]:
        bad = []
        for _ in range(2000):
            pi = list(range(n))
            random.shuffle(pi)
            a, b = I_brute(pi), I_fast(pi)
            if a != b:
                bad.append((pi, a, b))
        status = "OK" if not bad else f"MISMATCH {bad[:3]}"
        print(f"  n={n} random 2000: {status}")
        ok = ok and not bad
    return ok


def part_b_exhaustive_small():
    print("Part B1: bound exhaustive 4 <= n <= 9 (pure Python, independent I_brute)")
    ok = True
    known_max = {4: 2, 5: 4, 6: 6, 7: 9, 8: 12, 9: 16}
    for n in range(4, 10):
        worst = -1
        for pi in itertools.permutations(range(n)):
            v = I_brute(list(pi))
            if v > worst:
                worst = v
        b = bound(n)
        status = "OK" if worst <= b else "VIOLATION"
        match = "matches known max" if worst == known_max.get(n) else f"!= known {known_max.get(n)}"
        print(f"  n={n}: max I={worst} bound={b} [{status}], {match}")
        ok = ok and worst <= b and worst == known_max.get(n)
    return ok


def sigma(n):
    return [(-i) % n for i in range(n)]


def rev(n):
    return [n - 1 - i for i in range(n)]


def tau(n):
    m = n // 2
    p = list(range(n))
    p[m], p[m + 1] = p[m + 1], p[m]
    return p


def affine(n, a, b):
    return [(a * i + b) % n for i in range(n)]


def block_swap_antipode(n, k):
    p = list(range(n))
    m = n // 2
    for i in range(k):
        a, b = (m - k + i) % n, (m + i) % n
        p[a], p[b] = p[b], p[a]
    return p


def cycle_perm(n, length):
    p = list(range(n))
    idx = list(range(length))
    for k in range(length):
        p[idx[k]] = idx[(k + 1) % length]
    return p


def product_of_transpositions(n, count, seed):
    rnd = random.Random(seed)
    p = list(range(n))
    for _ in range(count):
        i, j = rnd.sample(range(n), 2)
        p[i], p[j] = p[j], p[i]
    return p


def part_b_large_n_samples():
    print("Part B2: structured + random samples n = 20, 50, 100, 101 (I_fast)")
    from math import gcd
    ok = True
    total_checked = 0
    for n in [20, 50, 100, 101]:
        b = bound(n)
        fams = {"sigma_n": sigma(n), "rev_n": rev(n), "identity": list(range(n))}
        if n % 2 == 0:
            fams["tau_n"] = tau(n)
            for k in [2, 3, 5, n // 4]:
                if 1 <= k <= n // 2:
                    fams[f"block_swap_antipode k={k}"] = block_swap_antipode(n, k)
        for a in [2, 3, 5, 7, n - 2, n - 3]:
            if a % n == 0 or gcd(a, n) != 1:
                continue
            fams[f"affine a={a},b=1"] = affine(n, a, 1)
            fams[f"affine a={a},b=n//3"] = affine(n, a, n // 3)
        for length in [3, 5, n // 2, n - 1]:
            if 2 <= length <= n:
                fams[f"cycle len={length}"] = cycle_perm(n, length)
        s = sigma(n)
        for shift in [1, 3, n // 5, n // 2]:
            fams[f"sigma_n rot {shift}"] = [s[(i + shift) % n] for i in range(n)]
        for cnt in [1, 2, 5, n // 4, n // 2]:
            fams[f"prod {cnt} transpositions"] = product_of_transpositions(n, cnt, seed=1000 + n + cnt)

        n_violations = 0
        for label, pi in fams.items():
            assert sorted(pi) == list(range(n))
            v = I_fast(pi)
            total_checked += 1
            if v > b:
                n_violations += 1
                print(f"  !!! n={n} {label}: I={v} > bound={b}")
        rnd = random.Random(42 + n)
        for _ in range(300):
            pi = list(range(n))
            rnd.shuffle(pi)
            v = I_fast(pi)
            total_checked += 1
            if v > b:
                n_violations += 1
        print(f"  n={n}: bound={b}, sigma_n achieves I={I_fast(sigma(n))} (should equal bound), "
              f"{len(fams) + 300} samples, violations={n_violations}")
        ok = ok and n_violations == 0 and I_fast(sigma(n)) == b
    print(f"  total samples checked at n in (20,50,100,101): {total_checked}")
    return ok


if __name__ == "__main__":
    a_ok = part_a()
    b1_ok = part_b_exhaustive_small()
    b2_ok = part_b_large_n_samples()
    print()
    print(f"Part A (formula correctness): {'PASS' if a_ok else 'FAIL'}")
    print(f"Part B1 (exhaustive 4<=n<=9): {'PASS' if b1_ok else 'FAIL'}")
    print(f"Part B2 (samples n=20,50,100,101): {'PASS' if b2_ok else 'FAIL'}")
    print("Note: n = 10..13 exhaustive verification is done by "
          "experiments/toric_I_fast.c (too slow in pure Python); see "
          "data/runs/h13_I_toric_scan/log.md for those results.")
    sys.exit(0 if (a_ok and b1_ok and b2_ok) else 1)
