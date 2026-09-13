"""h13_i_cut_restrictions-1.0 — H13-I proof-strategy screen.

For every pi (n small, exhaustive), compares I(pi) = min over all n^2 double
cuts (q, c) of the inversion count of the relabelled line (definition:
docs/notes/h13_line_model.md §0) against two restricted-search candidates
that would, if they worked, give a short proof of H13-I
(I(pi) <= floor((n-1)^2/4) for all pi, all n):

  * diagonal cuts only: c = q (n candidates instead of n^2);
  * position-only rotation: c fixed at a single canonical value, q varies
    (n candidates) — the classical "min inversions over cyclic rotations"
    quantity.

Usage: python3 h13_i_cut_restrictions.py [nmax]
"""
import itertools
import sys
import time


def bound(n):
    return ((n - 1) ** 2) // 4


def _rotated_line(pi, n, q, c):
    order_pos = [(i - q - 1) % n for i in range(n)]
    pos_at_rank = [0] * n
    for i in range(n):
        pos_at_rank[order_pos[i]] = i
    return [(pi[pos_at_rank[r]] - c - 1) % n for r in range(n)]


def _inversions(w):
    n = len(w)
    return sum(1 for a in range(n) for b in range(a + 1, n) if w[a] > w[b])


def I_full(pi, n):
    return min(
        _inversions(_rotated_line(pi, n, q, c))
        for q in range(n)
        for c in range(n)
    )


def I_diagonal(pi, n):
    return min(_inversions(_rotated_line(pi, n, q, q)) for q in range(n))


def I_position_only(pi, n, c=0):
    return min(_inversions(_rotated_line(pi, n, q, c)) for q in range(n))


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    for n in range(3, nmax + 1):
        t0 = time.time()
        B = bound(n)
        max_full = -1
        diag_fail = 0
        pos_fail = 0
        total = 0
        for pi in itertools.permutations(range(n)):
            total += 1
            pi = list(pi)
            full = I_full(pi, n)
            max_full = max(max_full, full)
            if full > B:
                print(f"COUNTEREXAMPLE to H13-I at n={n}: pi={pi}, I={full} > {B}")
            if I_diagonal(pi, n) > B:
                diag_fail += 1
            if I_position_only(pi, n) > B:
                pos_fail += 1
        dt = time.time() - t0
        print(
            f"n={n}: max I={max_full} bound={B} "
            f"diagonal_cuts_fail={diag_fail}/{total} "
            f"position_only_fail={pos_fail}/{total} "
            f"time={dt:.1f}s"
        )


if __name__ == "__main__":
    main()
