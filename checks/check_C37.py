"""Independent audit of C37: the double-cut inversion-grid gradient identity.

Claim (C37, CLAIMS.md): for a permutation pi of {0,...,n-1} and its inverse
invpi, define g(a, b) = inv(w) where w_j = (pi((a+j) mod n) - b) mod n for
j = 0..n-1 (the inversion count of the double cut (a, b), see
docs/notes/h13_line_model.md). Then:

    g(a+1, b) - g(a, b) = n - 1 - 2 * ((pi(a)      - b) mod n)   (mod n indices)
    g(a, b+1) - g(a, b) = n - 1 - 2 * ((invpi(b)   - a) mod n)

Proof sketch: g(a+1, b) is obtained from g(a, b) by moving the front element
of the b-shifted, a-rotated sequence (value m = (pi(a)-b) mod n) to the back;
this element loses its m inversions as the smallest-so-far element on the
left and gains n-1-m inversions as the new maximum on the right, a net
change of n-1-2m. The b-direction identity is the same argument applied to
the position where value b currently sits (dually, via invpi).

This script reimplements g(a, b) by direct brute force (independent of
experiments/toric_inversions_scan.c and toric_inversions_sample.c, which use
the identity to make exhaustive/large-n search feasible) and checks the
identity cell-by-cell for every permutation of every n in range.
"""
import itertools
import sys

N_MAX = 8  # exhaustive over all n! permutations up to this n


def inv_count(seq):
    n = len(seq)
    return sum(1 for i in range(n) for j in range(i + 1, n) if seq[i] > seq[j])


def brute_grid(pi, n):
    grid = [[0] * n for _ in range(n)]
    for a in range(n):
        v = pi[a:] + pi[:a]
        for b in range(n):
            w = [(x - b) % n for x in v]
            grid[a][b] = inv_count(w)
    return grid


def identity_grid(pi, n):
    invpi = [0] * n
    for i, v in enumerate(pi):
        invpi[v] = i
    grid = [[0] * n for _ in range(n)]
    grid[0][0] = inv_count(pi)
    for a in range(n - 1):
        grid[a + 1][0] = grid[a][0] + (n - 1 - 2 * pi[a])
    for a in range(n):
        for b in range(n - 1):
            jstar = (invpi[b] - a) % n
            grid[a][b + 1] = grid[a][b] + (n - 1 - 2 * jstar)
    return grid


def main():
    mismatches = 0
    checked = 0
    for n in range(2, N_MAX + 1):
        for pi in itertools.permutations(range(n)):
            pi = list(pi)
            checked += 1
            g1 = brute_grid(pi, n)
            g2 = identity_grid(pi, n)
            if g1 != g2:
                mismatches += 1
                if mismatches <= 5:
                    print(f"MISMATCH n={n} pi={pi}")
        print(f"n={n}: checked all {n}! permutations, mismatches so far: {mismatches}")

    print(f"TOTAL: {checked} permutations checked, {mismatches} mismatches")
    if mismatches:
        print("FAIL")
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
