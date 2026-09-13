"""Finite checks for C28 (docs/proofs/C28_reflections.md): H10 on reflections.

For pi(i) = h - i mod n every f_c is the reflection with parameter h' = h + c.
Lemma 1-3 of C28 give closed formulas for K_c (number of transpositions),
F_c, S_c by the parity of h' and the residue of n mod 4; Lemma 4 gives the
choice of c in {0, 1} and the value of the bound. This script compares the
formulas with the construction module (constructions/strict_upper.py,
strict_upper-1.0, core oracle-1.0) on ALL reflections and ALL shifts for
4 <= n <= NMAX, and checks the inequality of C28 with the code's R_c.

Each assert names the lemma it tests. The general argument is in the proof
text; these checks test the formulas and the implementation, not the proof.

Usage: python3 checks/check_C28.py [--nmax 60]
Output: data/runs/check_C28/report.json
"""

import argparse
import json
import os
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
for sub in ("oracle", "constructions", "experiments"):
    sys.path.insert(0, os.path.join(ROOT, sub))

from moves import CORE_VERSION  # noqa: E402
import strict_upper as su  # noqa: E402
from h10_scan import bound_R, B  # noqa: E402


def P(n):
    return n * n // 4


def formulas(n, hp):
    """Lemma 1-3: (K, F, S) of the reflection f(i) = hp - i."""
    if n % 2 == 1:
        K, F = (n - 1) // 2, P(n)
    elif hp % 2 == 0:
        K = (n - 2) // 2
        F = P(n) if n % 4 == 0 else P(n) - 1
    else:
        K = n // 2
        F = P(n) if n % 4 == 0 else P(n) + 1
    return K, F, 3 * K


def lemma4_choice(n, h):
    """Lemma 4: shift c in {0, 1} and the guaranteed value of 2F_c - S_c + H(c) - B_n."""
    if n % 2 == 1:
        return 1, 0
    if n % 4 == 0:
        return (0, 0) if h % 2 == 1 else (1, -1)
    # n = 2 mod 4
    return (1, 0) if h % 2 == 1 else (0, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=60)
    a = ap.parse_args()
    t0 = time.time()
    checked = 0
    plus_one_R = []
    for n in range(4, a.nmax + 1):
        Bn = B(n)
        for h in range(n):
            pi = tuple((h - i) % n for i in range(n))
            per_c = [bound_R(pi, c) for c in range(n)]
            for c, st in enumerate(per_c):
                hp = (h + c) % n
                K, F, S = formulas(n, hp)
                assert st["K"] == K, ("Lemma 1: K", n, h, c)
                assert st["F"] == F, ("Lemma 2: F", n, h, c)
                assert st["S"] == S, ("Lemma 3: S = 3K", n, h, c)
                assert st["R"] <= st["H"], ("R_c <= H(c)", n, h, c)
                checked += 1
            c0, val = lemma4_choice(n, h)
            st = per_c[c0]
            assert 2 * st["F"] - st["S"] + st["H"] - Bn == val, ("Lemma 4 value", n, h)
            assert val <= 1, ("Lemma 4: <= B_n + 1", n, h)
            # exactness of the H-route minimum (remark after Lemma 4)
            minH = min(2 * s["F"] - s["S"] + s["H"] for s in per_c) - Bn
            expected = 1 if (n % 4 == 2 and h % 2 == 0) else (0 if n % 2 == 1 or h % 2 == 1 else -1)
            assert minH == expected, ("H-route minimum", n, h, minH, expected)
            minR = min(s["bound"] for s in per_c) - Bn
            assert minR <= 1, ("C28 with R_c", n, h, minR)
            if minR == 1:
                plus_one_R.append((n, h))
    report = {
        "claim": "C28", "core": CORE_VERSION, "construction": su.CONSTRUCTION_VERSION,
        "range": f"4<=n<={a.nmax}, all reflections, all shifts",
        "checked_shifts": checked, "elapsed_s": round(time.time() - t0, 1),
        "R_route_plus_one_cases": plus_one_R,
        "result": "PASS",
    }
    out = os.path.join(ROOT, "data", "runs", "check_C28")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "report.json"), "w") as f:
        json.dump(report, f, indent=1)
    print(f"C28 PASS: {checked} (n, h, c) triples, 4<=n<={a.nmax}, {report['elapsed_s']} s")
    print("R-route minimum = B_n + 1 at (n, h):", plus_one_R)


if __name__ == "__main__":
    main()
