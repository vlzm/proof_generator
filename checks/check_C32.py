"""check_C32-1.0 (core oracle-1.0, construction strict_upper-1.0).

Verifies C32 (route bound via the largest carrier gap): for every pi in S_n,
c in Z_n, R_c <= 2(n - g_max(c)), where g_max(c) is the largest cyclic gap in
carriers(f_c) u {0,c} and R_c is the length of the carrier_route variant
(constructions/strict_upper.py:carrier_route). The bound is checked against
the actual route length returned by the audited carrier_route implementation,
exhaustively over all pi and all c for 4 <= n <= N (N given on the command
line).

Secondary report (not a proof, a measurement): coverage of the sharpened H10
sufficient condition (SC'), replacing R_c <= H(c) (C30/lemma 10) with
R_c <= 2(n - g_max) in the same inequality

    3 K_c + 2(P - F_c) >= floor(n/2) - 1 + 2(n - g_max(c))            (SC')

Existence of a c satisfying (SC') implies H10 for that pi (C29 gives
S_c >= 3 K_c, C32 gives R_c <= 2(n-g_max); substitute into lemma 8's
equivalence). (SC') is strictly sharper than (SC) of C30 (2(n-g_max) <= H(c)
always, see proof of C32), so its coverage is a lower bound on how much of
the general case these two lemmas alone can close; the gap that remains
CONJECTURED is reported per n, not swept under a sample.

Usage: python3 checks/check_C32.py --exhaustive N [--out DIR]
"""
import argparse
import itertools
import json
import os
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
sys.path.insert(0, os.path.join(ROOT, "constructions"))

from moves import CORE_VERSION  # noqa: E402
from strict_upper import shift_word, CONSTRUCTION_VERSION  # noqa: E402

VERSION = "check_C32-1.0"


def P_of(n):
    return n * n // 4


def gmax(pts, n):
    """Largest cyclic gap in a point set on Z_n; n if the set has <= 1 point."""
    pts = sorted(set(p % n for p in pts))
    m = len(pts)
    if m <= 1:
        return n
    best = 0
    for i in range(m):
        gap = (pts[(i + 1) % m] - pts[i]) % n
        if gap == 0:
            gap = n
        best = max(best, gap)
    return best


def check_n(n):
    P = P_of(n)
    half_floor = n // 2
    total_pairs = 0
    lemma_violations = []
    sc_hits = 0
    sc_prime_hits = 0
    n_perms = 0
    for pi in itertools.permutations(range(n)):
        n_perms += 1
        ok_sc = False
        ok_sc_prime = False
        for c in range(n):
            sw = shift_word(list(pi), c, route="carrier")
            st = sw["stats"]
            F_c, S_c, K_c = st["F_c"], st["S_c"], st["n_cycles"]
            R_c = st["route_len"]
            H_c = st["H"]
            pts = list(st["carriers"]) + [0, c % n]
            g = gmax(pts, n)
            bound = 2 * (n - g)
            total_pairs += 1
            if R_c > bound:
                lemma_violations.append(
                    {"pi": list(pi), "c": c, "R_c": R_c, "bound": bound, "gmax": g})
            lhs = 3 * K_c + 2 * (P - F_c)
            if lhs >= half_floor - 1 + H_c:
                ok_sc = True
            if lhs >= half_floor - 1 + bound:
                ok_sc_prime = True
        if ok_sc:
            sc_hits += 1
        if ok_sc_prime:
            sc_prime_hits += 1
    return {
        "n": n, "n_perms": n_perms, "pairs_checked": total_pairs,
        "lemma_violations": lemma_violations,
        "sc_coverage_pct": 100.0 * sc_hits / n_perms,
        "sc_prime_coverage_pct": 100.0 * sc_prime_hits / n_perms,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exhaustive", type=int, required=True,
                     help="check all pi for 4 <= n <= this value")
    ap.add_argument("--out", default=os.path.join(ROOT, "data", "runs", "check_C32"))
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    results = []
    t0 = time.time()
    for n in range(4, args.exhaustive + 1):
        tn0 = time.time()
        r = check_n(n)
        r["seconds"] = round(time.time() - tn0, 3)
        results.append(r)
        status = "PASS" if not r["lemma_violations"] else "FAIL"
        print(f"n={r['n']}: {r['n_perms']} pi, {r['pairs_checked']} (pi,c) pairs, "
              f"{r['seconds']}s, lemma {status}, "
              f"SC coverage {r['sc_coverage_pct']:.1f}%, "
              f"SC' coverage {r['sc_prime_coverage_pct']:.1f}%")

    report = {
        "version": VERSION, "core": CORE_VERSION, "construction": CONSTRUCTION_VERSION,
        "range": f"4<=n<={args.exhaustive}", "total_seconds": round(time.time() - t0, 3),
        "results": results,
    }
    with open(os.path.join(args.out, "report.json"), "w") as f:
        json.dump(report, f, indent=2)

    any_violation = any(r["lemma_violations"] for r in results)
    print()
    print("C32 (R_c <= 2(n-g_max)):", "FAIL" if any_violation else
          f"PASS, exhaustive 4<=n<={args.exhaustive}")
    if any_violation:
        sys.exit(1)


if __name__ == "__main__":
    main()
