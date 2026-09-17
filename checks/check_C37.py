"""check_C37 — independent exhaustive checker for the pair-corner formula
(docs/proofs/C37_pair_corner_formula.md).

For every permutation pi of {0,...,n-1} at 4 <= n <= 7, and every pair of
positions i1 < i2, verifies

    f(d,e) = n(d+e) - 2*d*e   (d = i2-i1, e = (pi(i2)-pi(i1)) mod n)

against a direct count, over all n^2 corners (a,b), of corners where the
pair is discordant. Independent implementation: re-derives "discordant"
from the raw definition (linear order after rotating positions by a and
values by b), not by importing experiments/h13_pair_corner.py.

PASS criterion: zero mismatches over all pairs of all permutations at
4 <= n <= 7 (24 + 120 + 720 + 5040 = 5904 permutations).
"""
import itertools
import time


def is_discordant(i1, i2, y1, y2, a, b, n):
    j1 = (i1 - a) % n
    j2 = (i2 - a) % n
    w1 = (y1 - b) % n
    w2 = (y2 - b) % n
    return (j1 < j2) != (w1 < w2)


def direct_count(i1, i2, y1, y2, n):
    return sum(
        1
        for a in range(n)
        for b in range(n)
        if is_discordant(i1, i2, y1, y2, a, b, n)
    )


def formula(d, e, n):
    return n * (d + e) - 2 * d * e


def main():
    t0 = time.time()
    total_perms = 0
    total_pairs = 0
    mismatches = 0
    for n in range(4, 8):
        for pi in itertools.permutations(range(n)):
            total_perms += 1
            for i1 in range(n):
                for i2 in range(i1 + 1, n):
                    d = i2 - i1
                    e = (pi[i2] - pi[i1]) % n
                    fd = formula(d, e, n)
                    direct = direct_count(i1, i2, pi[i1], pi[i2], n)
                    total_pairs += 1
                    if fd != direct:
                        mismatches += 1
                        print(f"MISMATCH n={n} pi={pi} i1={i1} i2={i2} "
                              f"d={d} e={e} formula={fd} direct={direct}")
        print(f"n={n}: {n*(n-1)//2} pairs x {__import__('math').factorial(n)} "
              f"perms checked, elapsed {time.time()-t0:.1f}s")
    print(f"TOTAL permutations: {total_perms}, pairs: {total_pairs}, "
          f"mismatches: {mismatches}")
    print("PASS" if mismatches == 0 else "FAIL")


if __name__ == "__main__":
    main()
