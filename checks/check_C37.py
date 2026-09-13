"""Checker for C37 (docs/proofs/C37_toric_correlation.md).

Independent re-implementation: inversions are recomputed from the H13-I
definition of the line (w_j = (pi(q+1+j) - (q+1-c)) mod n) by direct pair
counting; sig, S, E are rebuilt here and not imported from the experiment.

  Part A (Theorem A, identity).  For every pi and every cut (q, c),
      4 n^2 inv(q,c) = n^2 (n^2 - 1) - E(pi) - 2 n S(q, (q-c) mod n).
    Exhaustive over all pi and all n^2 cuts for 2 <= n <= AMAX (default 6),
    plus random pi for n up to RMAX (default 12).
    Also checks the intermediate lemmas 1-3 of the proof:
      sig(q-i) = 2 x_i - (n-1),  S = 4 sum_j j w_j - n(n-1)^2,
      inv + (2/n) sum_j j w_j constant over the n^2 cuts.
  Part B (Theorem B, energy bound).  E(pi) >= 2 n^2 (n-1) - n^2 (n^2-1)/3 with
    equality exactly on the n reflections pi_h(i) = h - i: exhaustive for
    2 <= n <= BMAX (default 8), random for larger n, reflections up to 200.
  Part C (Corollary C).  avg_{q,c} inv = (n^2-1)/4 - E/(4n^2) <= (n-1)(2n-1)/6
    and I(pi) <= floor((n-1)(2n-1)/6); the exhaustive maxima max_pi I(pi) are
    printed against floor((n-1)^2/4) (H13-I target) for 4 <= n <= BMAX.
  Part D (remarks 3-5 of the proof).  The equivalent form of H13-I
    (2 n max S + E >= 2 n thr), the reflection values max S and E, and the
    closed form S = C((q+s-h) mod n) for reflections.

Usage: python3 checks/check_C37.py [--amax 6] [--bmax 8] [--rmax 12]
Output: data/runs/check_C37/report.json, report.md; exit code 0 iff PASS.
"""

import argparse
import itertools
import json
import os
import random
import sys
import time
from fractions import Fraction

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUTDIR = os.path.join(ROOT, "data", "runs", "check_C37")


def sig(n, t):
    return n - 1 - 2 * (t % n)


def line(perm, q, c):
    n = len(perm)
    return [(perm[(q + 1 + j) % n] - (q + 1 - c)) % n for j in range(n)]


def inv_count(w):
    n = len(w)
    return sum(1 for j in range(n) for k in range(j + 1, n) if w[j] > w[k])


def S_at(perm, q, s):
    n = len(perm)
    return sum(sig(n, q - i) * sig(n, s - perm[i]) for i in range(n))


def E_of(perm):
    n = len(perm)
    return sum(sig(n, j - i) * sig(n, perm[j] - perm[i])
               for i in range(n) for j in range(n))


def E_bound(n):
    return 2 * n * n * (n - 1) - n * n * (n * n - 1) // 3


def reflections(n):
    return [tuple((h - i) % n for i in range(n)) for h in range(n)]


def check_perm_identity(perm):
    """Theorem A + lemmas 1-3 for one permutation.  Returns list of errors."""
    n = len(perm)
    errs = []
    E = E_of(perm)
    phi = None
    for q in range(n):
        for c in range(n):
            s = (q - c) % n
            w = line(perm, q, c)
            inv = inv_count(w)
            Sv = S_at(perm, q, s)
            # lemma 1
            for i in range(n):
                x = (i - q - 1) % n
                y = (perm[i] - (q + 1 - c)) % n
                if sig(n, q - i) != 2 * x - (n - 1):
                    errs.append(("L1x", perm, q, c, i))
                if sig(n, s - perm[i]) != 2 * y - (n - 1):
                    errs.append(("L1y", perm, q, c, i))
            # lemma 2
            mom = sum(j * w[j] for j in range(n))
            if Sv != 4 * mom - n * (n - 1) ** 2:
                errs.append(("L2", perm, q, c))
            # lemma 3
            cur = Fraction(inv) + Fraction(2 * mom, n)
            if phi is None:
                phi = cur
            elif cur != phi:
                errs.append(("L3", perm, q, c))
            # theorem A
            if 4 * n * n * inv != n * n * (n * n - 1) - E - 2 * n * Sv:
                errs.append(("A", perm, q, c))
    return errs


def part_A(amax, rmax, rng, trials):
    res = {"exhaustive": {}, "random": {}, "errors": []}
    for n in range(2, amax + 1):
        errs = []
        for perm in itertools.permutations(range(n)):
            errs += check_perm_identity(list(perm))
        res["exhaustive"][n] = {"perms": len(list(itertools.permutations(range(n)))),
                                "errors": len(errs)}
        res["errors"] += errs[:3]
    for n in range(amax + 1, rmax + 1):
        errs = []
        for _ in range(trials):
            p = list(range(n))
            rng.shuffle(p)
            errs += check_perm_identity(p)
        res["random"][n] = {"trials": trials, "errors": len(errs)}
        res["errors"] += errs[:3]
    return res


def part_B(bmax, rmax, rng, trials):
    res = {"exhaustive": {}, "random": {}, "reflections": {}}
    for n in range(2, bmax + 1):
        bnd = E_bound(n)
        eq = set()
        below = 0
        for perm in itertools.permutations(range(n)):
            E = E_of(perm)
            if E < bnd:
                below += 1
            elif E == bnd:
                eq.add(perm)
        res["exhaustive"][n] = {"bound": bnd, "below": below,
                                "equality_is_reflections": eq == set(reflections(n))}
    for n in range(bmax + 1, rmax + 1):
        bnd = E_bound(n)
        below = 0
        eq_non_refl = 0
        refl = set(reflections(n))
        for _ in range(trials):
            p = list(range(n))
            rng.shuffle(p)
            E = E_of(p)
            if E < bnd:
                below += 1
            if E == bnd and tuple(p) not in refl:
                eq_non_refl += 1
        res["random"][n] = {"trials": trials, "below": below,
                            "equality_off_reflections": eq_non_refl}
    for n in list(range(4, 41)) + [60, 100, 200]:
        bnd = E_bound(n)
        ok = all(E_of(list(p)) == bnd for p in reflections(n))
        res["reflections"][n] = ok
    return res


def part_C(bmax):
    res = {}
    for n in range(4, bmax + 1):
        bnd = ((n - 1) * (2 * n - 1)) // 6
        target = (n - 1) ** 2 // 4
        max_I = 0
        viol_avg = 0
        viol_I = 0
        for perm in itertools.permutations(range(n)):
            E = E_of(perm)
            smax = max(S_at(perm, q, s) for q in range(n) for s in range(n))
            num = n * n * (n * n - 1) - E - 2 * n * smax
            if num % (4 * n * n) != 0:
                viol_I += 1000
            I = num // (4 * n * n)
            # I from brute force too
            I2 = min(inv_count(line(perm, q, c))
                     for q in range(n) for c in range(n))
            if I != I2:
                viol_I += 1000
            max_I = max(max_I, I)
            if I > bnd:
                viol_I += 1
            avg = Fraction(sum(inv_count(line(perm, q, c))
                               for q in range(n) for c in range(n)), n * n)
            if avg != Fraction(n * n - 1, 4) - Fraction(E, 4 * n * n):
                viol_avg += 1000
            if avg > Fraction((n - 1) * (2 * n - 1), 6):
                viol_avg += 1
        res[n] = {"proved_bound": bnd, "max_I": max_I, "h13i_target": target,
                  "violations_I": viol_I, "violations_avg": viol_avg}
    return res


def part_D(bmax):
    """Remarks: equivalent form of H13-I and the reflection closed form."""
    res = {"equivalence_failures": 0, "reflection_formula_failures": 0,
           "reflection_maxS": {}, "checked": 0}
    for n in range(4, bmax + 1):
        thr = n * (n // 2) + n * (n - 1) // 2
        peak = -n * (n * n - 1) // 3 + 2 * n * (n * n // 4)
        for perm in itertools.permutations(range(n)):
            E = E_of(perm)
            smax = max(S_at(perm, q, s) for q in range(n) for s in range(n))
            I = min(inv_count(line(perm, q, c))
                    for q in range(n) for c in range(n))
            lhs = (I <= (n - 1) ** 2 // 4)
            rhs = (2 * n * smax + E >= 2 * n * thr)
            if lhs != rhs:
                res["equivalence_failures"] += 1
            res["checked"] += 1
    for n in range(4, 21):
        peak = -n * (n * n - 1) // 3 + 2 * n * (n * n // 4)
        res["reflection_maxS"][n] = peak
        for h in range(n):
            p = [(h - i) % n for i in range(n)]
            smax = max(S_at(p, q, s) for q in range(n) for s in range(n))
            if smax != peak:
                res["reflection_formula_failures"] += 1
            for q in range(n):
                for s in range(n):
                    m = (q + s - h) % n
                    C = (-n * (n * n - 1) // 3
                         + 2 * n * (m + 1) * (n - 1 - m))
                    if S_at(p, q, s) != C:
                        res["reflection_formula_failures"] += 1
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amax", type=int, default=6)
    ap.add_argument("--bmax", type=int, default=7)
    ap.add_argument("--rmax", type=int, default=12)
    ap.add_argument("--trials", type=int, default=20)
    ap.add_argument("--seed", type=int, default=37)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    t0 = time.time()

    rep = {"args": vars(args)}
    rep["A"] = part_A(args.amax, args.rmax, rng, args.trials)
    rep["B"] = part_B(args.bmax, args.rmax, rng, args.trials)
    rep["C"] = part_C(args.bmax)
    rep["D"] = part_D(min(args.bmax, 7))
    rep["seconds"] = round(time.time() - t0, 1)

    ok = True
    ok &= all(d["errors"] == 0 for d in rep["A"]["exhaustive"].values())
    ok &= all(d["errors"] == 0 for d in rep["A"]["random"].values())
    ok &= all(d["below"] == 0 and d["equality_is_reflections"]
              for d in rep["B"]["exhaustive"].values())
    ok &= all(d["below"] == 0 and d["equality_off_reflections"] == 0
              for d in rep["B"]["random"].values())
    ok &= all(rep["B"]["reflections"].values())
    ok &= all(d["violations_I"] == 0 and d["violations_avg"] == 0
              for d in rep["C"].values())
    ok &= rep["D"]["equivalence_failures"] == 0
    ok &= rep["D"]["reflection_formula_failures"] == 0
    rep["status"] = "PASS" if ok else "FAIL"

    os.makedirs(OUTDIR, exist_ok=True)
    with open(os.path.join(OUTDIR, "report.json"), "w") as f:
        json.dump(rep, f, indent=1, sort_keys=True)
    lines = [f"# check_C37 — {rep['status']} ({rep['seconds']} s)", "",
             "Part A (Theorem A + lemmas 1-3): exhaustive "
             + ", ".join(f"n={n}: {d['errors']} errors"
                         for n, d in sorted(rep["A"]["exhaustive"].items()))
             + "; random " + ", ".join(f"n={n}: {d['errors']} errors"
                                       for n, d in sorted(rep["A"]["random"].items())),
             "",
             "Part B (Theorem B): " + ", ".join(
                 f"n={n}: below {d['below']}, equality=reflections "
                 f"{d['equality_is_reflections']}"
                 for n, d in sorted(rep["B"]["exhaustive"].items())),
             "  reflections hit the bound for n up to 200: "
             + str(all(rep["B"]["reflections"].values())), "",
             "Part C (Corollary C):"]
    for n, d in sorted(rep["C"].items()):
        lines.append(f"  n={n}: proved bound {d['proved_bound']}, max_pi I "
                     f"{d['max_I']}, H13-I target {d['h13i_target']}, "
                     f"violations {d['violations_I']}/{d['violations_avg']}")
    lines += ["", f"Part D: equivalence failures {rep['D']['equivalence_failures']}"
                  f" on {rep['D']['checked']} permutations; reflection closed-form "
                  f"failures {rep['D']['reflection_formula_failures']}"]
    with open(os.path.join(OUTDIR, "report.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
