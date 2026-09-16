"""H13-I averaging probe (session 9): does "best q, then pigeonhole over c"
prove I(pi) <= floor((n-1)^2/4) for all pi?

For a fixed position-cut q, sum_{j<k} e_jk(q) = n * (average over c of the
number of non-inverted pairs).  By pigeonhole max_c concordant(q, c) >=
ceil(sum_jk e_jk(q) / n).  This script computes, for every pi, the best such
guarantee over all q and compares it to the true I(pi) (found by exhaustive
n^2 search) and to the target floor(n^2/4) non-inversions.

Result (see docs/notes/h13_line_model.md, session 9, section 7): the
guarantee falls short of the target on pi = id by n - 3, growing with n,
while the true I(id) = 0 (trivial: q = n-1, c = 0 sorts it).  So the
averaging technique itself is too weak to prove H13-I, independent of
whether the theorem is true; this is a defect of the method, not evidence
against the conjecture.  Version h13_averaging_probe-1.0.
"""

import itertools


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def exact_I(pi):
    """I(pi) = min over all n^2 cuts (q, c) of the number of inversions."""
    n = len(pi)
    best = n * n
    for q in range(n):
        for c in range(n):
            shift = (q + 1 - c) % n
            w = [(pi[(q + 1 + j) % n] - shift) % n for j in range(n)]
            best = min(best, inversions(w))
    return best


def sum_e_jk(pi, q):
    """sum_{j<k} e_jk(q) = n * avg_c concordant(q, c); see docstring above."""
    n = len(pi)
    tau = [pi[(q + 1 + j) % n] for j in range(n)]
    return sum((tau[k] - tau[j]) % n for j in range(n) for k in range(j + 1, n))


def best_q_pigeonhole_guarantee(pi):
    """max over q of the pigeonhole guarantee on max_c concordant(q, c)."""
    n = len(pi)
    best = 0
    for q in range(n):
        total = sum_e_jk(pi, q)
        best = max(best, -(-total // n))  # ceil(total / n)
    return best


def main(nmax=8):
    for n in range(4, nmax + 1):
        worst_shortfall = None
        worst_pi = None
        for pi in itertools.permutations(range(n)):
            target = n * n // 4
            guarantee = best_q_pigeonhole_guarantee(list(pi))
            shortfall = target - guarantee
            if worst_shortfall is None or shortfall > worst_shortfall:
                worst_shortfall = shortfall
                worst_pi = pi
        print(
            f"n={n} target(non-inv)={n * n // 4} "
            f"worst shortfall={worst_shortfall} pi={worst_pi} "
            f"true I(pi)={exact_I(list(worst_pi))}"
        )


if __name__ == "__main__":
    main()
