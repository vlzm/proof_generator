"""check_C1-1.0 (core oracle-1.0, construction strict_upper-1.0).

Independent audit of C1 (`docs/incoming/lrx_2n3_proof_checked.md`, "An upper
bound for the LRX diameter"): D_n <= floor(A_n) at n >= 4, where
A_n = 2P + n - (2P+11n-11)/(3n), P = floor(n^2/4); consequence
D_n <= B_n + floor(4n/3) - 3.

The per-cycle local word construction of C1 sections 1-2 (shortest signed
paths, carrier choice at a max-load rise, contraction, disjoint XX
cancellations) is the same construction later reused and refined by N1
(`constructions/strict_upper.py`, already audited as strict_upper-1.0): the
word for shift c and the identity len == 2*F_c - S_c + H(c) are already
independently checked by `checks/check_C28.py`/`check_C32.py`. This checker
targets what is specific to C1 and not already covered by the N1 audit:

1. Inequality (5), `3 n S(C) >= 5(F(C)+k)`, for every nontrivial cycle C -
   this is C1's own per-cycle saving bound (weaker and structurally
   different from N1's C29, `S(C) >= 3`); checked exhaustively on every
   cycle produced by every (pi, c) at 4 <= n <= N.
2. The exact-arithmetic rounding chain: floor(R(n)) (the real bound of the
   averaging argument) equals the source code's own integer `stated_bound`,
   and is <= the theorem's second, weaker closed form
   B_n + floor(4n/3) - 3, checked for 4 <= n <= N_SCALAR by exact integer
   arithmetic (no floats).
3. Section 6's sharpness family C_r (n = 4r): k=8, t=4, E=0, F=3n-8, M=3,
   S=5, with equality in (5).
4. Cross-execution: the word constructed by `docs/incoming/lrx_2n3_check.py`
   (the source's own reference implementation) is replayed with
   `oracle/moves.py` (an independently audited engine, not the source
   file's own `apply_word`) and confirmed to sort pi within the stated
   integer bound, for every (pi, c) at 4 <= n <= N.

Usage: python3 checks/check_C1.py --exhaustive N [--scalar N_SCALAR] [--sharp R] [--out DIR]
"""
import argparse
import importlib.util
import itertools
import json
import os
import sys
import time
from fractions import Fraction

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
sys.path.insert(0, os.path.join(ROOT, "constructions"))

from moves import apply_word, identity, CORE_VERSION  # noqa: E402
from strict_upper import (cycles_of, cycle_data, shift_word,  # noqa: E402
                           CONSTRUCTION_VERSION, signed_step)

VERSION = "check_C1-1.0"

_spec = importlib.util.spec_from_file_location(
    "lrx_2n3_check", os.path.join(ROOT, "docs", "incoming", "lrx_2n3_check.py"))
_src = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = _src
_spec.loader.exec_module(_src)


def P_of(n):
    return n * n // 4


def B_of(n):
    return n * (n - 1) // 2


def R_exact(n):
    """Exact value of the real bound A_n (theorem, first display)."""
    P = P_of(n)
    return 2 * P + n - Fraction(2 * P + 11 * n - 11, 3 * n)


def floor_frac(x):
    return x.numerator // x.denominator


# ------------------------------------------------------------- 1. inequality (5)

def check_ineq5(n_max, out):
    total = 0
    violations = []
    tight = []
    for n in range(4, n_max + 1):
        for pi in itertools.permutations(range(n)):
            for c in range(n):
                f = [(pi[i] + c) % n for i in range(n)]
                for cyc in cycles_of(f):
                    cd = cycle_data(cyc, n)
                    F, k, S = cd["F"], cd["k"], cd["S"]
                    total += 1
                    lhs = 3 * n * S
                    rhs = 5 * (F + k)
                    if lhs < rhs:
                        violations.append({"n": n, "cycle": cyc, "F": F, "k": k,
                                            "M": cd["M"], "E": cd["E"], "S": S})
                    elif lhs == rhs:
                        tight.append((n, k, cd["M"], cd["E"]))
    report = {"cycles_checked": total, "violations": violations[:20],
              "n_violations": len(violations),
              "tight_cases_by_n": sorted(set(tight))}
    print(f"(5) 3nS(C)>=5(F+k): {total} cycle instances, 4<=n<={n_max}, "
          f"{len(violations)} violations, {len(set(tight))} distinct tight "
          f"(n,k,M,E) patterns")
    return report


# ------------------------------------------------------ 2. rounding / scalar chain

def check_rounding(n_scalar):
    mismatches = []
    unsound = []
    for n in range(4, n_scalar + 1):
        R = R_exact(n)
        fR = floor_frac(R)
        sb = _src.stated_bound(n)
        if fR != sb:
            mismatches.append({"n": n, "floor_R": fR, "stated_bound": sb})
        second = B_of(n) + (4 * n) // 3 - 3
        if fR > second:
            unsound.append({"n": n, "floor_R": fR, "second_bound": second})
    print(f"rounding chain: 4<=n<={n_scalar}, floor(R(n))==stated_bound(n) "
          f"mismatches={len(mismatches)}; floor(R(n))<=B_n+floor(4n/3)-3 "
          f"violations={len(unsound)}")
    return {"n_scalar": n_scalar, "mismatches": mismatches[:20],
            "unsound": unsound[:20]}


# --------------------------------------------------------------- 3. sharpness (S6)

def check_sharpness(r_max):
    results = []
    ok = True
    for r in range(2, r_max + 1):
        n = 4 * r
        cyc = [0, 2 * r - 1, r, 3 * r - 1, 2 * r, 4 * r - 1, 3 * r, r - 1]
        assert len(set(cyc)) == 8, (r, cyc)
        cd = cycle_data(cyc, n)
        exp_F = 3 * n - 8
        good = (cd["k"] == 8 and cd["t"] == 4 and cd["E"] == 0 and
                cd["F"] == exp_F and cd["M"] == 3 and cd["S"] == 5 and
                3 * n * cd["S"] == 5 * (cd["F"] + cd["k"]))
        ok = ok and good
        results.append({"r": r, "n": n, **cd, "matches_S6": good})
    print(f"section 6 sharpness family C_r: r=2..{r_max}, all match: {ok}")
    return {"ok": ok, "sample": results[:5]}


# ---------------------------------------------------- 4. cross-execution of words

def check_cross_execution(n_max):
    total = 0
    failures = []
    for n in range(4, n_max + 1):
        bound = _src.stated_bound(n)
        for pi in itertools.permutations(range(n)):
            best = None
            for c in range(n):
                word = _src.word_for_shift(list(pi), c)
                # replay with the independently audited oracle engine, not the
                # source file's own apply_word
                final = apply_word(tuple(pi), "".join(word))
                total += 1
                if final != identity(n):
                    failures.append({"n": n, "pi": pi, "c": c, "issue": "oracle replay failed"})
                    continue
                if best is None or len(word) < best:
                    best = len(word)
            if best is None or best > bound:
                failures.append({"n": n, "pi": pi, "issue": "bound exceeded",
                                  "best": best, "bound": bound})
    print(f"cross-execution: {total} (pi,c) words replayed with oracle/moves.py, "
          f"4<=n<={n_max}, {len(failures)} failures")
    return {"total": total, "failures": failures[:20], "n_failures": len(failures)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exhaustive", type=int, default=7,
                     help="exhaustive n range for (1) and (4): 4<=n<=this")
    ap.add_argument("--scalar", type=int, default=100000,
                     help="range for the exact rounding chain (2)")
    ap.add_argument("--sharp", type=int, default=25,
                     help="max r for the section-6 sharpness family (3), n=4r")
    ap.add_argument("--out", default=os.path.join(ROOT, "data", "runs", "check_C1"))
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    report = {
        "version": VERSION, "core": CORE_VERSION, "construction": CONSTRUCTION_VERSION,
        "ineq5": check_ineq5(args.exhaustive, args.out),
        "rounding": check_rounding(args.scalar),
        "sharpness": check_sharpness(args.sharp),
        "cross_execution": check_cross_execution(args.exhaustive),
    }
    report["seconds"] = round(time.time() - t0, 3)
    with open(os.path.join(args.out, "report.json"), "w") as f:
        json.dump(report, f, indent=2, default=str)

    fail = (report["ineq5"]["n_violations"] or report["rounding"]["mismatches"] or
            report["rounding"]["unsound"] or not report["sharpness"]["ok"] or
            report["cross_execution"]["n_failures"])
    print()
    print("check_C1:", "FAIL" if fail else f"PASS, exhaustive 4<=n<={args.exhaustive}, "
          f"scalar 4<=n<={args.scalar}, sharp r<={args.sharp}, {report['seconds']}s")
    if fail:
        sys.exit(1)


if __name__ == "__main__":
    main()
