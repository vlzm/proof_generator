"""Independent check of the closed-form proof docs/proofs/H10_reflections.md
(H10 / C27 restricted to reflections pi(i) = h - i mod n).

The proof claims exact closed-form values of F_c, S_c, H(c) for reflections at
c in {0, 1} (odd n needs only c = 1), an exact formula for the resulting
excess 2F_c - S_c + H(c) - B_n depending only on (n mod 4, parity of h), and
the corollary min_c[2F_c - S_c + R_c] <= B_n + 1 (route R_c <= H(c) always).

This script recomputes F_c, S_c, H(c) with constructions/strict_upper.py
(the trusted, independently-tested implementation of N1 sections 1-3) for a
range of n and ALL h, and checks them against the closed-form formulas of the
proof, then checks the derived excess bound, then separately re-derives the
same bound using the actual carrier-route construction (route="carrier"),
executing and verifying every produced word. Passing finite checks do not
replace the algebraic proof (AGENTS.md rule 5); the algebraic argument is
n-independent (the derivation eliminates n from the excess formula
algebraically) and is checked by full case coverage of the four
(n mod 4, h mod 2) classes plus the separate odd-n argument, not by taking
n large.

Usage: python3 checks/check_C27_reflections.py [--max-n-exhaustive 60]
       [--max-n-formula 400] [--max-n-carrier 40]
Output: data/runs/check_C27_reflections/report.json, summary.md
"""

import argparse
import json
import os
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
for sub in ("oracle", "constructions", "bounds"):
    sys.path.insert(0, os.path.join(ROOT, sub))

from moves import apply_word, freely_reduce, identity, CORE_VERSION  # noqa: E402
import known  # noqa: E402
import strict_upper as su  # noqa: E402

OUT_DIR = os.path.join(ROOT, "data", "runs", "check_C27_reflections")


class CheckFailure(Exception):
    pass


def _req(cond, msg, info=None):
    if not cond:
        raise CheckFailure(f"{msg}: {info!r}")


def reflection(n, h):
    return tuple((h - i) % n for i in range(n))


def raw_stats(pi, c):
    """F_c, S_c, H(c) recomputed directly from N1 sections 1-3 (cycle_data),
    independent of the aggregate stats dict in shift_word (cross-check)."""
    n = len(pi)
    f = [(pi[i] + c) % n for i in range(n)]
    cycles = su.cycles_of(f)
    F = sum(su.cycle_data(cyc, n)["F"] for cyc in cycles)
    S = sum(su.cycle_data(cyc, n)["S"] for cyc in cycles)
    return F, S, su.H_of(c, n), len(cycles)


# --------------------------------------------------------------------- Lemma 1 (odd n)

def closed_form_odd(n):
    """Proof section 2: F_c = P, S_c = 3(n-1)/2 for ALL c, when n is odd."""
    P = known.P(n)
    return P, 3 * (n - 1) // 2


def check_odd(max_n, log):
    """Lemma 1 (F_c, S_c constant for odd n, every c) and Theorem (excess = 0
    at c = 1) checked for every odd n in range and every h, every c."""
    checked_full, checked_c1 = 0, 0
    for n in range(5, max_n + 1, 2):
        Bn = known.target_diameter(n)
        F_expected, S_expected = closed_form_odd(n)
        for h in range(n):
            pi = reflection(n, h)
            for c in range(n):
                F, S, H, K = raw_stats(pi, c)
                _req(F == F_expected, "odd n: F_c != P for all c", (n, h, c, F, F_expected))
                _req(S == S_expected, "odd n: S_c != 3(n-1)/2 for all c", (n, h, c, S, S_expected))
                checked_full += 1
            # theorem: c = 1 alone already gives excess 0 (exact B_n)
            F, S, H, K = raw_stats(pi, 1)
            _req(H == n - 1, "odd n: H(1) != n-1", (n, H))
            excess = 2 * F - S + H - Bn
            _req(excess == 0, "odd n: 2F_1-S_1+H(1) != B_n", (n, h, excess))
            checked_c1 += 1
    log(f"odd n 5..{max_n}: F_c=P and S_c=3(n-1)/2 at every (h,c) -- {checked_full} triples; "
        f"2F_1-S_1+H(1)=B_n exactly at every h -- {checked_c1} cases")
    return {"range": f"5<=n<={max_n} step 2 (all odd)", "triples_full": checked_full, "cases_c1": checked_c1}


# --------------------------------------------------------------------- Lemma 2-4 (even n)

def closed_form_even(n, h, c):
    """Proof sections 3-4: exact F_c, S_c, H(c) for even n at c in {0,1},
    as a function of n mod 4 and the parities of h, k = h + c."""
    _req(n % 2 == 0 and c in (0, 1), "closed_form_even domain", (n, c))
    P = known.P(n)
    k_even = (h + c) % 2 == 0
    S = 3 * (n - 2) // 2 if k_even else 3 * n // 2
    if n % 4 == 0:
        F = P
    else:  # n % 4 == 2
        F = P - 1 if k_even else P + 1
    H = n if c == 0 else n - 1
    return F, S, H


def check_even(max_n, log):
    checked = 0
    excess_classes = {}   # (n%4, h%2) -> set of excess-at-{0,1} pairs, expect singleton
    for n in range(4, max_n + 1, 2):
        Bn = known.target_diameter(n)
        for h in range(n):
            pi = reflection(n, h)
            best = None
            excs = []
            for c in (0, 1):
                F, S, H, K = raw_stats(pi, c)
                Fe, Se, He = closed_form_even(n, h, c)
                _req((F, S, H) == (Fe, Se, He), "even n closed form mismatch",
                     (n, h, c, (F, S, H), (Fe, Se, He)))
                exc = 2 * F - S + H - Bn
                excs.append(exc)
                best = exc if best is None else min(best, exc)
                checked += 1
            key = (n % 4, h % 2)
            excess_classes.setdefault(key, set()).add(tuple(excs))
            _req(best <= 1, "even n: min_{c in {0,1}} excess > 1", (n, h, best))
    for key, vals in sorted(excess_classes.items()):
        _req(len(vals) == 1, "excess pair not constant within (n%4,h%2) class", (key, vals))
    log(f"even n 4..{max_n}: closed-form F_c,S_c,H(c) match construction at c=0,1 for every h "
        f"-- {checked} (h,c) pairs; excess pairs by class: "
        + ", ".join(f"n%4={k[0]},h%2={k[1]}:{sorted(v)[0]}" for k, v in sorted(excess_classes.items())))
    return {"range": f"4<=n<={max_n} step 2 (all even)",
            "pairs_checked": checked,
            "excess_by_class": {f"n%4={k[0]},h%2={k[1]}": sorted(v)[0] for k, v in sorted(excess_classes.items())}}


# --------------------------------------------------------------------- corollary: actual words, both routes

def check_words(max_n, log, route):
    """Build and execute the real word (both routes) for every reflection at
    every n in range, confirm it sorts pi and len - B_n matches the bound."""
    worst = -10 ** 9
    worst_ex = None
    n_words = 0
    for n in range(4, max_n + 1):
        Bn = known.target_diameter(n)
        for h in range(n):
            pi = reflection(n, h)
            best_len = None
            for c in range(n):
                sw = su.shift_word(pi, c, route)
                _req(apply_word(pi, sw["word"]) == identity(n), "word does not sort pi", (pi, c, route))
                best_len = sw["stats"]["len"] if best_len is None else min(best_len, sw["stats"]["len"])
                n_words += 1
            exc = best_len - Bn
            _req(exc <= 1, f"H10 route={route}: min_c length exceeds B_n+1", (n, h, exc))
            if exc > worst:
                worst = exc
                worst_ex = (n, h)
    log(f"route={route}: all reflections 4<=n<={max_n}, all shifts, words executed and verified to sort "
        f"({n_words} words); worst min_c(len)-B_n = {worst} at {worst_ex}")
    return {"range": f"4<=n<={max_n}", "route": route, "words_checked": n_words,
            "worst_excess": worst, "worst_example_n_h": worst_ex}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n-exhaustive", type=int, default=60,
                     help="exhaustive closed-form check range (all h, c in {0,1} even / all c odd)")
    ap.add_argument("--max-n-formula", type=int, default=400,
                     help="reserved for future larger-range spot checks (unused: proof is n-independent)")
    ap.add_argument("--max-n-carrier", type=int, default=40,
                     help="range for executing actual carrier-route words (route='carrier')")
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)
    t0 = time.time()
    logs = []

    def log(msg):
        print(msg, flush=True)
        logs.append(msg)

    report = {"core": CORE_VERSION, "construction": su.CONSTRUCTION_VERSION,
              "claim": "H10/C27 restricted to reflections pi(i)=h-i mod n "
                       "(docs/proofs/H10_reflections.md)",
              "command": " ".join(sys.argv)}
    report["odd_n"] = check_odd(args.max_n_exhaustive, log)
    report["even_n"] = check_even(args.max_n_exhaustive, log)
    report["words_n1_route"] = check_words(args.max_n_carrier, log, "n1")
    report["words_carrier_route"] = check_words(args.max_n_carrier, log, "carrier")
    report["seconds"] = round(time.time() - t0, 1)
    report["log"] = logs

    with open(os.path.join(OUT_DIR, "report.json"), "w") as f:
        json.dump(report, f, indent=1)
    with open(os.path.join(OUT_DIR, "summary.md"), "w") as f:
        f.write("# check_C27_reflections — summary\n\n")
        for line in logs:
            f.write(f"- {line}\n")
        f.write(f"\nTotal time: {report['seconds']} s.\n")
    print("PASS", report["seconds"], "s")


if __name__ == "__main__":
    main()
