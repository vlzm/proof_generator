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

from moves import CORE_VERSION, freely_reduce, apply_word, identity  # noqa: E402
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


def check_part2(nmax):
    """C28 part 2 (docs/proofs/C28_reflections.md, 'Часть 2'): after free
    reduction, the carrier-route word has length <= B_n on every reflection.
    Only the exceptional case of part 1 (n = 2 mod 4, h even) needs checking;
    cases 1-4 already give <= B_n before reduction (Lemma 4 value <= 0)."""
    worst = None
    pairs = 0
    for n in range(4, nmax + 1):
        if n % 4 != 2:
            continue
        Bn = B(n)
        for h in range(0, n, 2):
            pi = tuple((h - i) % n for i in range(n))
            best = None
            for c in range(n):
                sw = su.shift_word(pi, c, route="carrier")
                red = freely_reduce(sw["word"])
                assert apply_word(pi, red) == identity(n), ("part 2 word must sort pi", n, h, c)
                if best is None or len(red) < best:
                    best = len(red)
            pairs += 1
            excess = best - Bn
            assert excess <= 0, ("C28 part 2: <= B_n after reduction", n, h, excess)
            if worst is None or excess > worst[0]:
                worst = (excess, n, h)
    return {"pairs_checked": pairs, "worst_excess": worst[0] if worst else None,
            "worst_at": worst[1:] if worst else None}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=60)
    ap.add_argument("--part2-nmax", type=int, default=40,
                     help="range for the (slower) part-2 reduction check")
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
    t1 = time.time()
    part2 = check_part2(a.part2_nmax)
    t2 = time.time()
    report = {
        "claim": "C28", "core": CORE_VERSION, "construction": su.CONSTRUCTION_VERSION,
        "range": f"4<=n<={a.nmax}, all reflections, all shifts",
        "checked_shifts": checked, "elapsed_s": round(t1 - t0, 1),
        "R_route_plus_one_cases": plus_one_R,
        "part2": {"range": f"4<=n<={a.part2_nmax}, n=2 mod 4, all even h",
                  "elapsed_s": round(t2 - t1, 1), **part2},
        "result": "PASS",
    }
    out = os.path.join(ROOT, "data", "runs", "check_C28")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "report.json"), "w") as f:
        json.dump(report, f, indent=1)
    print(f"C28 part 1 PASS: {checked} (n, h, c) triples, 4<=n<={a.nmax}, {report['elapsed_s']} s")
    print("R-route minimum = B_n + 1 at (n, h):", plus_one_R)
    print(f"C28 part 2 PASS: {part2['pairs_checked']} (n, h) pairs, "
          f"4<=n<={a.part2_nmax}, worst excess {part2['worst_excess']}, {report['part2']['elapsed_s']} s")


if __name__ == "__main__":
    main()
