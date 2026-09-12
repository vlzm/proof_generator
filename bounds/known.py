"""Known bounds and their exact integer evaluation (PLAN.md §9.10).

Statuses and dependencies live in CLAIMS.md; a function existing here does
not certify its claim. All floors use exact integer arithmetic (no float).
"""

from fractions import Fraction
from math import isqrt


def target_diameter(n):
    """B_n = n(n-1)/2 — the conjectured diameter (C3)."""
    return n * (n - 1) // 2


def P(n):
    return n * n // 4


def A_n(n):
    """Rational A_n from the baseline construction (C1)."""
    p = P(n)
    return 2 * p + n - Fraction(2 * p + 11 * n - 11, 3 * n)


def baseline_bound_2n3(n):
    """floor(A_n): claimed D_n <= floor(A_n) for n >= 4 (C1, CLAIMED)."""
    a = A_n(n)
    return a.numerator // a.denominator


def gap_bound(n):
    """Gap-file improvements (C4/C10/C11, CLAIMED)."""
    b = baseline_bound_2n3(n) - 1
    if n >= 138:
        b -= 1
    if n >= 852:
        b -= 1
    return b


def ceil_sqrt(x):
    """Exact ceiling of sqrt for non-negative integers."""
    r = isqrt(x)
    return r if r * r == x else r + 1


def strict_upper_bound(n):
    """U_n = floor(R_n) from N1 (C16, CLAIMED).

    R_n = (a - b*sqrt(z)) / d with
    z = 3n^2 - 4n + 1 - 4P, a = 36nP + 25n^2 - 43n - 12P + 36,
    b = 7n, d = 18n; hence U_n = (a - ceil_sqrt(b^2 z)) // d.
    floor((a - b*sqrt(z))/d) = (a - ceil(b*sqrt(z))) // d holds because for
    integer a, d > 0 and real y: floor((a - y)/d) = (a - ceil(y)) // d, and
    ceil(b*sqrt(z)) = ceil_sqrt(b^2 z) (exact when b^2 z is a perfect square,
    since then b*sqrt(z) is an integer; otherwise b*sqrt(z) is irrational and
    both ceilings agree).
    """
    p = P(n)
    z = 3 * n * n - 4 * n + 1 - 4 * p
    a = 36 * n * p + 25 * n * n - 43 * n - 12 * p + 36
    b = 7 * n
    d = 18 * n
    return (a - ceil_sqrt(b * b * z)) // d


def strict_upper_bound_simple(n):
    """C17: B_n + floor(kappa*n) - 2 (odd n), -1 (even n); kappa=(31-7*sqrt2)/18.

    floor(kappa*n) = floor((31n - 7n*sqrt(2)) / 18)
                   = (31n - ceil_sqrt(98 n^2)) // 18   (98 n^2 never a square).
    """
    fk = (31 * n - ceil_sqrt(98 * n * n)) // 18
    return target_diameter(n) + fk - (2 if n % 2 else 1)


CLAIMED_BOUNDS = (
    ("C1", 4, baseline_bound_2n3),
    ("C4/C10/C11", 4, gap_bound),
    ("C16", 4, strict_upper_bound),
    ("C17", 4, strict_upper_bound_simple),
)  # C1/C4/C10/C11 are CLAIMED; C16/C17 are PROVED (see PROVED_BOUNDS)

# Locally VERIFIED exact diameters (certified tables in data/tables).
VERIFIED_EXACT = {4: 6, 5: 10, 6: 15, 7: 21, 8: 28, 9: 36, 10: 45, 11: 55,
                  12: 66}


def best_claimed(n):
    vals = [(f(n), cid) for cid, n0, f in CLAIMED_BOUNDS if n >= n0]
    v, cid = min(vals)
    return v, cid


# Bounds PROVED by local audit (CLAIMS.md): C16 (U_n) and C17, both n >= 4,
# audit of 2026-09-12 in docs/notes/strict_proof_audit.md.
PROVED_BOUNDS = (
    ("C16", 4, strict_upper_bound),
    ("C17", 4, strict_upper_bound_simple),
)


def certified_bound(n):
    """Confirmed bound: min of VERIFIED exact data (finite range) and
    locally PROVED general bounds (C16/C17). C1/C4/C10/C11 stay CLAIMED
    and are excluded here."""
    vals = [(f(n), cid) for cid, n0, f in PROVED_BOUNDS if n >= n0]
    if n in VERIFIED_EXACT:
        vals.append((VERIFIED_EXACT[n], "VERIFIED (certified BFS table)"))
    if not vals:
        return None, "no locally certified bound for this n"
    return min(vals)


if __name__ == "__main__":
    # control values from PLAN.md §9.10
    assert [strict_upper_bound(n) for n in (4, 7, 20, 100)] == [9, 27, 211, 5065]
    assert [baseline_bound_2n3(n) for n in (4, 5, 6, 7, 8)] == [8, 13, 19, 26, 35]
    print("control values OK")
    print(f"{'n':>4} {'B_n':>6} {'C1':>6} {'gap':>6} {'U_n':>6} {'C17':>6} "
          f"{'best(claimed)':>14} {'certified':>10}")
    for n in (4, 5, 6, 7, 8, 9, 10, 20, 50, 100, 137, 138, 851, 852, 1000):
        v, cid = best_claimed(n)
        cert, _ = certified_bound(n)
        print(f"{n:>4} {target_diameter(n):>6} {baseline_bound_2n3(n):>6} "
              f"{gap_bound(n):>6} {strict_upper_bound(n):>6} "
              f"{strict_upper_bound_simple(n):>6} {v:>8} ({cid:>4}) "
              f"{str(cert):>10}")
