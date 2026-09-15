"""toric_I_recurrence_check.py — session 9, H13-I.

Brute-force verification, for small n, of two exact identities used by
experiments/toric_I_fast.c to fill the n x n grid inv(s, t) (s = position
rotation, t = value rotation; inv(s, t) is the inversion count of the line
w_j = (pi[(s+j) mod n] - t) mod n, j = 0..n-1 — the same operational
definition of the double cut as experiments/line_profile.c) in O(n^2) instead
of O(n^4) time per permutation:

  inv(s+1, t) - inv(s, t) = (n-1) - 2*((pi[s]    - t) mod n)     (s-recurrence)
  inv(s, t+1) - inv(s, t) = (n-1) - 2*((piinv[t] - s) mod n)     (t-recurrence)

and the exact pair-counting formula for the sum over all n^2 cuts:

  sum_{s,t} inv(s,t) = sum over unordered pairs {p,q} of points (i, pi(i))
                        of f(a,b),  f(a,b) = n(a+b) - 2ab
                        = n^2/2 - 2(a-n/2)(b-n/2),
  where a = (x2-x1) mod n, b = (y2-y1) mod n for the two points (x1,y1),(x2,y2)
  of the pair (either orientation of the pair gives the same f, since
  f(n-a,n-b) = f(a,b)).

See docs/notes/h13_line_model.md section 7 for the derivation and discussion
(these recurrences are the same content as C33's Lemma C, applied k=1 step at
a time and now proved as exact identities, not just necessary conditions at a
local optimum).

Usage: python3 experiments/toric_I_recurrence_check.py --nmax 7
Version toric_I_recurrence_check-1.0.
"""

import argparse
import itertools


def inv_count(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def line(pi, n, s, t):
    return [(pi[(s + j) % n] - t) % n for j in range(n)]


def brute_grid(pi, n):
    return {(s, t): inv_count(line(pi, n, s, t)) for s in range(n) for t in range(n)}


def check_recurrences(pi, n, grid):
    for s in range(n):
        for t in range(n):
            w0 = (pi[s % n] - t) % n
            lhs = grid[((s + 1) % n, t)] - grid[(s, t)]
            rhs = (n - 1) - 2 * w0
            assert lhs == rhs, ("s-recurrence failed", pi, s, t, lhs, rhs)
    piinv = [0] * n
    for i in range(n):
        piinv[pi[i]] = i
    for s in range(n):
        for t in range(n):
            k = (piinv[t % n] - s) % n
            lhs = grid[(s, (t + 1) % n)] - grid[(s, t)]
            rhs = (n - 1) - 2 * k
            assert lhs == rhs, ("t-recurrence failed", pi, s, t, lhs, rhs)


def check_pair_formula(pi, n, grid):
    total = sum(grid.values())
    pts = [(i, pi[i]) for i in range(n)]
    s2 = 0
    for i in range(n):
        for j in range(i + 1, n):
            x1, y1 = pts[i]
            x2, y2 = pts[j]
            a = (x2 - x1) % n
            b = (y2 - y1) % n
            s2 += n * (a + b) - 2 * a * b
    assert total == s2, ("pair formula failed", pi, total, s2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=2)
    ap.add_argument("--nmax", type=int, default=7)
    args = ap.parse_args()
    count = 0
    for n in range(args.nmin, args.nmax + 1):
        for pi in itertools.permutations(range(n)):
            pi = list(pi)
            grid = brute_grid(pi, n)
            check_recurrences(pi, n, grid)
            check_pair_formula(pi, n, grid)
            count += 1
    print("all checks passed (s-recurrence, t-recurrence, pair formula f(a,b)); "
          "permutations tested:", count, "for n in [%d, %d]" % (args.nmin, args.nmax))


if __name__ == "__main__":
    main()
