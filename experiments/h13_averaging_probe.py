"""H13-I probe (session 9): test whether a one-dimensional averaging argument
over double cuts (q, c) suffices to prove I(pi) <= floor((n-1)^2/4).

For a permutation pi of Z_n and a position cut a (0 <= a < n), let
s_j = pi((a+j) mod n) for j = 0..n-1 (line with position order fixed by a).
A(a) = average over value shifts b of inv(a,b), where inv(a,b) counts
inversions of the sequence ((s_j - b) mod n)_j.

We check max_pi min_a A(a) against floor((n-1)^2/4): if min_a A(a) stayed
within the bound, picking the best a and then any b at or below the average
(pigeonhole) would prove H13-I without an explicit choice of b. Exhaustive
over all pi, all n; see docs/notes/h13_line_model.md ss7.1 for the verdict.
"""
import itertools
import sys


def inversions(seq):
    n = len(seq)
    c = 0
    for i in range(n):
        for j in range(i + 1, n):
            if seq[i] > seq[j]:
                c += 1
    return c


def A_of_a(pi, a, n):
    s = [pi[(a + j) % n] for j in range(n)]
    total = 0
    for j in range(n):
        for k in range(j + 1, n):
            total += (s[k] - s[j]) % n
    return total / n


def min_over_ab_bruteforce(pi, n):
    best = None
    for a in range(n):
        s = [pi[(a + j) % n] for j in range(n)]
        for b in range(n):
            t = [(x - b) % n for x in s]
            iv = inversions(t)
            if best is None or iv < best:
                best = iv
    return best


def main():
    nmin = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    print(f"{'n':>2} {'bound':>6} {'max min_a A(a)':>15} {'argmax pi':>0}")
    for n in range(nmin, nmax + 1):
        bound = (n - 1) ** 2 // 4
        worst_minA = -1.0
        worst_pi = None
        worst_control = None
        for pi in itertools.permutations(range(n)):
            minA = min(A_of_a(list(pi), a, n) for a in range(n))
            if minA > worst_minA:
                worst_minA = minA
                worst_pi = pi
                worst_control = min_over_ab_bruteforce(list(pi), n)
        print(f"{n:>2} {bound:>6} {worst_minA:>15.4f} {worst_pi} "
              f"(control: true I(pi) = {worst_control})")


if __name__ == "__main__":
    main()
