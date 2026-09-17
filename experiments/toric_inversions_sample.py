"""toric_inversions_sample.py — structured/random samples of I(pi) at large n
(AGENTS.md rule 9 ladder step: n = 20, 50, 100, 101 and neighbourhoods).

I(pi) = min over double cuts (q, c) of the number of inversions of the
relabelled line (definition: docs/notes/h13_line_model.md §0). Computed with
the same O(n^2) row/column recursion as experiments/toric_inversions_fast.c
(SAMPLED here, not exhaustive: this checks specific pi, not all of S_n).

Families: both reflections (h = 0, 1), the antipodal transposition (n even),
invertible affine maps pi(i) = a*i+b mod n, uniformly random permutations
(fixed seed), identity, rev_n. This is error-finding, not a proof (rule 9).

Usage: python3 experiments/toric_inversions_sample.py
Version toric_inversions_sample-1.0.
"""

import math
import random


def I_fast(pi):
    n = len(pi)
    pinv = [0] * n
    for i, v in enumerate(pi):
        pinv[v] = i

    def inv0(seq):
        m = len(seq)
        c = 0
        for i in range(m):
            for j in range(i + 1, m):
                if seq[i] > seq[j]:
                    c += 1
        return c

    base = inv0(pi)
    row0 = [0] * n
    row0[0] = base
    for b in range(n - 1):
        r = pinv[b]
        row0[b + 1] = row0[b] + (n - 1) - 2 * r
    best = 10 ** 9
    for b in range(n):
        cur = row0[b]
        if cur < best:
            best = cur
        for a in range(n - 1):
            r = pi[a] - b
            if r < 0:
                r += n
            cur = cur + (n - 1) - 2 * r
            if cur < best:
                best = cur
    return best


def check(name, pi, n, rows):
    bound = ((n - 1) ** 2) // 4
    val = I_fast(pi)
    status = "OK" if val <= bound else "VIOLATION"
    rows.append((n, name, val, bound, status))
    print(f"n={n:4d} {name:30s} I={val:6d} bound={bound:6d} {status}")


def main():
    rows = []
    for n in (20, 50, 100, 101):
        check("reflection h=0", [(-i) % n for i in range(n)], n, rows)
        check("reflection h=1", [(1 - i) % n for i in range(n)], n, rows)
        if n % 2 == 0:
            pi = list(range(n))
            m = n // 2
            pi[m], pi[m + 1] = pi[m + 1], pi[m]
            check("antipodal transposition", pi, n, rows)
        for a, b in [(2, 0), (3, 1), (5, 7), (n - 2, 3)]:
            if math.gcd(a, n) != 1:
                continue
            check(f"affine a={a},b={b}", [(a * i + b) % n for i in range(n)], n, rows)
        random.seed(1)
        for trial in range(3):
            pi = list(range(n))
            random.shuffle(pi)
            check(f"random#{trial}", pi, n, rows)
        check("identity", list(range(n)), n, rows)
        check("rev_n", [n - 1 - i for i in range(n)], n, rows)
    violations = [r for r in rows if r[4] == "VIOLATION"]
    print(f"\ntotal checks: {len(rows)}, violations: {len(violations)}")


if __name__ == "__main__":
    main()
