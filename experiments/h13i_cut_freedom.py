"""h13i_cut_freedom-1.0 -- two negative results about candidate proof routes
for H13-I (I(pi) <= floor((n-1)^2/4) for all pi, all n), checked exhaustively
for small n. Both confirm and quantify remarks already on file in
docs/notes/h13_line_model.md section 6 ("averaging does not work"); this
script is the reproducible source for the numbers.

1. average_concordant(pi): the average, over all n^2 double cuts (q, c), of
   the number of concordant (non-inverted) pairs of the relabelled line.
   If this average were >= floor(n^2/4) for every pi, H13-I would follow
   immediately (some cut is at least as good as the average). It is not:
   the average is well *below* the threshold precisely on the permutations
   that are extremal for I(pi) (the reflections), so plain averaging over
   all n^2 cuts cannot prove H13-I in either direction.

2. best_value_rotation_only(pi): fix the position order (q = n-1, i.e. no
   position rotation) and take the best of only the n value-rotations c.
   This is a strictly weaker search (n choices instead of n^2) and is
   insufficient starting at n = 7: it exceeds floor((n-1)^2/4) by 1, showing
   the position-rotation freedom is not redundant, i.e. H13-I genuinely
   needs the double (position, value) freedom, not just one of the two.

Usage: python3 h13i_cut_freedom.py --nmax 8
"""
import argparse
import itertools


def inv_count(w):
    n = len(w)
    c = 0
    for i in range(n):
        for j in range(i + 1, n):
            if w[i] > w[j]:
                c += 1
    return c


def bound_I(n):
    return ((n - 1) ** 2) // 4


def total_inv_over_all_cuts(pi):
    n = len(pi)
    total = 0
    for q in range(n):
        for c in range(n):
            w = [(pi[(q + 1 + j) % n] - (q + 1 - c)) % n for j in range(n)]
            total += inv_count(w)
    return total


def best_value_rotation_only(pi):
    n = len(pi)
    best = None
    for c in range(n):
        w = [(pi[i] - c) % n for i in range(n)]
        iv = inv_count(w)
        if best is None or iv < best:
            best = iv
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=4)
    ap.add_argument("--nmax", type=int, default=8)
    args = ap.parse_args()

    print("n  min_avg_concordant  threshold_floor(n^2/4)  averaging_suffices")
    for n in range(args.nmin, args.nmax + 1):
        threshold = (n * n) // 4
        worst_avg = None
        worst_pi = None
        for pi in itertools.permutations(range(n)):
            total = total_inv_over_all_cuts(list(pi))
            avg_conc = n * (n - 1) / 2 - total / (n * n)
            if worst_avg is None or avg_conc < worst_avg:
                worst_avg = avg_conc
                worst_pi = pi
        print(f"{n}  {worst_avg:.4f}  {threshold}  {worst_avg >= threshold}  argmin={worst_pi}")

    print()
    print("n  max_over_pi(best_value_rotation_only)  bound_floor((n-1)^2/4)  sufficient")
    for n in range(args.nmin, args.nmax + 1):
        b = bound_I(n)
        worst = 0
        worst_pi = None
        for pi in itertools.permutations(range(n)):
            v = best_value_rotation_only(list(pi))
            if v > worst:
                worst = v
                worst_pi = pi
        print(f"{n}  {worst}  {b}  {worst <= b}  argmax={worst_pi}")


if __name__ == "__main__":
    main()
