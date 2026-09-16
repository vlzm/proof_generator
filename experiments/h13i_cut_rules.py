"""h13i_cut_rules.py — H13-I candidate cut-selection rules (h13i_cut_rules-1.0).

H13-I (docs/notes/h13_line_model.md §6): every pi has a double cut (q, c) with
inv(w) <= floor((n-1)^2/4). This script tests whether two simple *single-shot*
rules for choosing q (independent of c, then best c for that q) already
achieve the bound, as a step toward a constructive proof. Both fail (off by
exactly 1) already at n = 7, 8 — negative result, recorded in
data/runs/h13i_cut_rules/report.md and docs/notes/h13_line_model.md.

w_j = pi[(q+1+j) % n] - (q+1-c) mod n, j = 0..n-1 (same convention as
experiments/line_profile.c).
"""
import itertools
import sys


def inv_count(w):
    n = len(w)
    c = 0
    for i in range(n):
        for j in range(i + 1, n):
            if w[i] > w[j]:
                c += 1
    return c


def line_inv(pi, q, c, n):
    shift = (q + 1 - c) % n
    w = [(pi[(q + 1 + j) % n] - shift) % n for j in range(n)]
    return inv_count(w)


def target(n):
    return ((n - 1) ** 2) // 4


def worst_case(n, rule, limit=None):
    """rule(pi, n) -> chosen q. Returns (worst achieved value, worst pi)."""
    worst, worst_pi, cnt = -1, None, 0
    for pi in itertools.permutations(range(n)):
        cnt += 1
        if limit and cnt > limit:
            break
        q = rule(pi, n)
        best_c_val = min(line_inv(pi, q, c, n) for c in range(n))
        if best_c_val > worst:
            worst, worst_pi = best_c_val, pi
    return worst, worst_pi


def rule_min_inv0(pi, n):
    """Choose q minimizing inv at fixed shift c = 0."""
    return min(range(n), key=lambda q: line_inv(pi, q, 0, n))


def rule_min_S(pi, n):
    """Choose q minimizing S(q) = sum_c inv(q, c) (average over c)."""
    return min(range(n), key=lambda q: sum(line_inv(pi, q, c, n) for c in range(n)))


RULES = {"min_inv0": rule_min_inv0, "min_S": rule_min_S}

if __name__ == "__main__":
    n = int(sys.argv[1])
    limit = int(sys.argv[2]) if len(sys.argv) > 2 else None
    for name, rule in RULES.items():
        worst, worst_pi = worst_case(n, rule, limit)
        tgt = target(n)
        print(f"rule={name} n={n} target={tgt} worst_achieved={worst} "
              f"pi={worst_pi} PASS={worst <= tgt}")
