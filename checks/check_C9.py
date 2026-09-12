"""Check of C9 via the N1 construction (docs/proofs/C9_sigma_upper.md).

Claim: for every n >= 4 the N1 word (constructions/strict_upper.py) for
sigma_n with shift c* = 0 (n = 0 mod 4) or c* = 1 (otherwise) has length
exactly B_n = n(n-1)/2, hence d(sigma_n) <= B_n.

Finite part: for 4 <= n <= N the word is built, executed with the reference
moves and its length, N_X and N_rot are recorded; the closed forms of
Lemma 2 (F_{c*}, S_{c*}, H(c*)) are asserted. The general statement is the
algebraic proof in docs/proofs/C9_sigma_upper.md; this script checks the
implementation and the closed forms, it does not replace the proof.
"""

import argparse
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "oracle"))
sys.path.insert(0, os.path.join(HERE, "..", "constructions"))
from moves import apply_word, sigma, CORE_VERSION  # noqa: E402
import strict_upper as SU  # noqa: E402


def c_star(n):
    return 0 if n % 4 == 0 else 1


def closed_forms(n):
    """(F, S, H, number of transpositions) of Lemma 2 for c = c_star(n)."""
    P = n * n // 4
    if n % 2 == 1:
        m = (n - 1) // 2
        return P, 3 * m, n - 1, m
    m = n // 2
    if m % 2 == 0:
        return P, 3 * m, n, m
    return P - 1, 3 * (m - 1), n - 1, m - 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=40)
    ap.add_argument("--report")
    args = ap.parse_args()
    rows = []
    t0 = time.time()
    for n in range(4, args.max_n + 1):
        pi = sigma(n)
        c = c_star(n)
        s = SU.shift_word(pi, c)
        Bn = n * (n - 1) // 2
        F, S, H, ntr = closed_forms(n)
        assert s["F_c"] == F and s["S_c"] == S and s["H"] == H, "Lemma 2"
        assert len(s["cycles"]) == ntr and all(len(C) == 2 for C in s["cycles"])
        for lw in s["locals"]:
            st = lw["stats"]
            assert st["M"] == 2 and st["E"] == 0 and st["S"] == 3, \
                "Lemma 1: every transposition is non-antipodal, S = 3"
        assert s["len"] == 2 * F - S + H == Bn, "Lemma 3: length B_n"
        assert apply_word(pi, s["word"]) == tuple(range(n)), "word sorts"
        rows.append({"n": n, "c": c, "len": s["len"], "B_n": Bn,
                     "N_X": s["N_X"], "N_rot": s["N_rot"],
                     "reduced_len": s["reduced_len"],
                     "floor((n-1)^2/4)": (n - 1) ** 2 // 4,
                     "floor(n^2/4)": n * n // 4})
        print(f"n={n:3d} c*={c} len={s['len']:5d} = B_n  N_X={s['N_X']:5d} "
              f"(floor((n-1)^2/4)={(n-1)**2//4})  N_rot={s['N_rot']:5d} "
              f"(floor(n^2/4)={n*n//4})  reduced={s['reduced_len']}")
    print(f"C9 words verified for 4<=n<={args.max_n} in {time.time()-t0:.1f}s")
    if args.report:
        with open(args.report, "w") as f:
            json.dump({"core_version": CORE_VERSION,
                       "construction_version": SU.CONSTRUCTION_VERSION,
                       "command": " ".join(sys.argv),
                       "range": f"4<=n<={args.max_n}", "rows": rows},
                      f, indent=1)
    print("PASS")


if __name__ == "__main__":
    main()
