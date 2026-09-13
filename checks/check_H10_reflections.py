"""Independent checker for docs/proofs/H10_reflections.md (C28).

Two things are verified separately, both in exact integer arithmetic:

1. The closed-form formulas for F_c, S_c of a reflection pi(i) = h - i mod n
   (section 2 of the proof) against a from-scratch reimplementation of the
   construction's cycle statistics (not importing constructions/strict_upper.py,
   per AGENTS.md rule 3 for independent audits: re-derive from the definitions
   in PROBLEM.md section 5.2, not from the audited module).
2. The final inequality min_c [2F_c - S_c + R_c] <= B_n + 1 via the explicit
   strategy of section 3 (c in {0, 1}, R_0 <= n, R_1 <= H(1) = n-1), in exact
   integer arithmetic, for a large range of n (scale of the C16v scalar check).

Usage: python3 checks/check_H10_reflections.py
"""

import sys


def P(n):
    return n * n // 4


def target_diameter(n):
    return n * (n - 1) // 2


# ---------------------------------------------------------------- part 1: F_c, S_c
# Independent reimplementation of the cycle statistics of f(i) = (h - i) mod n,
# following PROBLEM.md section 5.2 definitions directly (signed shortest step,
# antipodal steps positive; F(C) = sum |d|, S(C) = 2M(C) + E(C) - 1).

def signed_step(a, b, n):
    d = (b - a) % n
    if 2 * d > n:
        d -= n
    return d


def reflection_stats(hp, n):
    """F_c, S_c, theta_c for f(i) = (hp - i) mod n, by brute force over all
    nontrivial 2-cycles (reflections have no cycles longer than 2)."""
    seen = [False] * n
    F = 0
    S = 0
    theta = 0
    fixed = 0
    for i in range(n):
        if seen[i]:
            continue
        j = (hp - i) % n
        if j == i:
            fixed += 1
            seen[i] = True
            continue
        assert not seen[j], (hp, n, i, j)
        seen[i] = seen[j] = True
        d1 = signed_step(i, j, n)
        d2 = signed_step(j, i, n)
        F += abs(d1) + abs(d2)
        antipodal = (n % 2 == 0 and d1 == n // 2 and d2 == n // 2)
        if antipodal:
            theta += 1
            M, E = 1, 2
        else:
            assert d1 == -d2, (hp, n, i, j, d1, d2)
            M, E = 2, 0
        S += 2 * M + E - 1
    return F, S, theta, fixed


def closed_form(n, hp):
    hp_parity = hp % 2
    if n % 2 == 1:
        return P(n), 3 * (n - 1) // 2
    n4 = n % 4
    if n4 == 0:
        F = P(n)
        S = (3 * n) // 2 - 3 if hp_parity == 0 else (3 * n) // 2
    else:
        F = P(n) - 1 if hp_parity == 0 else P(n) + 1
        S = (3 * n) // 2 - 3 if hp_parity == 0 else (3 * n) // 2
    return F, S


def check_formulas():
    n_full = list(range(4, 60)) + [99, 100, 101, 102, 199, 200, 201, 202]
    checked = 0
    for n in n_full:
        hp_range = range(n) if n <= 60 else range(4)  # representatives by parity for large n
        for hp in hp_range:
            F, S, theta, fixed = reflection_stats(hp, n)
            F_pred, S_pred = closed_form(n, hp)
            assert F == F_pred, ("F mismatch", n, hp, F, F_pred)
            assert S == S_pred, ("S mismatch", n, hp, S, S_pred)
            # cross-check theta/fixed against section 2.1 case table
            if n % 2 == 1:
                assert fixed == 1 and theta == 0
            else:
                m = n // 2
                if hp % 2 == 0:
                    assert fixed == 2
                    assert theta == (1 if m % 2 == 0 else 0)
                else:
                    assert fixed == 0
                    assert theta == (0 if m % 2 == 0 else 1)
            checked += 1
    print(f"part 1 OK: {checked} (n, h') pairs, formulas + fixed/theta case table match")


# ---------------------------------------------------------------- part 2: H10 bound
# Uses only the two literal values from N1 section 3 (H(0) = n, H(1) = n - 1),
# not a general H(c) formula, since the strategy only ever needs c in {0, 1}.

def check_bound(n_max):
    worst = None
    violations = 0
    for n in range(4, n_max + 1):
        Bn = target_diameter(n)
        if n % 2 == 1:
            F, S = closed_form(n, 1)  # any h' works, constant
            L = 2 * F - S
            bound = L + n  # c = 0, R_0 <= H(0) = n
            excess = bound - (Bn + 1)
        else:
            n4 = n % 4
            good_parity = 1 if n4 == 0 else 0
            excess = None
            for h_parity in (0, 1):
                if h_parity == good_parity:
                    c, Rbound = 0, n
                else:
                    c, Rbound = 1, n - 1
                hprime_parity = (h_parity + c) % 2
                F, S = closed_form(n, hprime_parity)
                L = 2 * F - S
                bound = L + Rbound
                e = bound - (Bn + 1)
                excess = e if excess is None else max(excess, e)
        if excess > 0:
            violations += 1
            print("VIOLATION at n =", n, "excess =", excess)
        if worst is None or excess > worst[0]:
            worst = (excess, n)
    print(f"part 2 OK: n = 4..{n_max}, worst (min_c bound) - (B_n+1) = {worst[0]} at n={worst[1]}, "
          f"violations = {violations}")
    assert violations == 0


if __name__ == "__main__":
    check_formulas()
    n_max = 100000
    if len(sys.argv) > 1:
        n_max = int(sys.argv[1])
    check_bound(n_max)
    print("H10 on reflections (C28): PASS")
