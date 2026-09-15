"""h13i_avgq_pairwise_check.py — exhaustive check of the closed form for
sum_q inv(v^{(q)}) used in data/runs/h13i_avgq/report.md (the "position-only"
half of the H13-I-avg decomposition; does not cover the harder min_c S_c
half). Version h13i_avgq_pairwise_check-1.0.

Claim: for pi a permutation of {0,...,n-1} and v^{(q)}_j = pi[(q+1+j) mod n],
    sum_{q=0}^{n-1} inv(v^{(q)}) = sum_{0<=i1<i2<=n-1} f(i1, i2),
    f(i1, i2) = (n - d) if pi[i1] > pi[i2] else d,   d = i2 - i1.

Usage: python3 h13i_avgq_pairwise_check.py [nmax]  (default nmax=8, exhaustive)
"""
import itertools
import sys


def inv(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def sum_q_inv(pi):
    n = len(pi)
    return sum(inv([pi[(q + 1 + j) % n] for j in range(n)]) for q in range(n))


def formula(pi):
    n = len(pi)
    s = 0
    for i1 in range(n):
        for i2 in range(i1 + 1, n):
            d = i2 - i1
            s += (n - d) if pi[i1] > pi[i2] else d
    return s


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    for n in range(2, nmax + 1):
        checked = 0
        for pi in itertools.permutations(range(n)):
            a, b = sum_q_inv(list(pi)), formula(list(pi))
            if a != b:
                print(f"MISMATCH n={n} pi={pi} sum_q_inv={a} formula={b}")
                return 1
            checked += 1
        print(f"n={n}: {checked} permutations, formula matches sum_q inv(v^(q)) exactly")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
