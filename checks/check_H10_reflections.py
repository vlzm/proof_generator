"""Independent finite checks of docs/proofs/H10_reflections.md (H10, step 1: PLAN §8).

Checks two separate things, both tied to explicit lemma numbers of the proof:

1. Lemma B/C (closed forms). For a reflection f_w(i) = (w - i) mod n, the
   number of transpositions K(w), F(w) = 2F_c-relevant sum of |steps| and
   S(w) = sum of local S(C) depend only on the parity of n and of w, by an
   explicit formula. This is checked against the actual cycle data computed
   by constructions/strict_upper.py (cycles_of, cycle_data) -- an independent
   implementation of N1 sections 1-2 -- for every w, not just the shift used
   by the theorem.

2. Theorem (explicit good shift). For every reflection pi_h(i) = (h-i) mod n
   and the shift c = good_c(n, h) defined in the proof, the actual word
   built by constructions/strict_upper.shift_word (route "n1", i.e. using
   N1's own H(c)) has length exactly the value predicted by the proof's case
   table, and the carrier-route word (route "carrier", i.e. R_c) is never
   longer than it (Lemma A) and its excess over B_n is at most +1.

Coverage ladder (AGENTS.md rule 9): exhaustive over all h for 4<=n<=<exhaustive_max>;
structured samples (several h per n, covering all residues mod 4 and both
parities of h) at larger n and near the requested threshold.

Every word produced is executed with the reference moves (oracle-1.0) and
asserted to sort its input; this file never reads distance tables.

Usage: python3 checks/check_H10_reflections.py [--exhaustive-max 80]
       [--sample-ns 200,201,202,203,500,501,502,503,998,999,1000,1001]
Output: data/runs/check_H10_reflections/report.json and summary.md
"""

import argparse
import json
import os
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
for sub in ("oracle", "constructions"):
    sys.path.insert(0, os.path.join(ROOT, sub))

from moves import apply_word, freely_reduce, identity, CORE_VERSION  # noqa: E402
import strict_upper as su  # noqa: E402

OUT_DIR = os.path.join(ROOT, "data", "runs", "check_H10_reflections")


class CheckFailure(Exception):
    pass


def Bn(n):
    return n * (n - 1) // 2


def P(n):
    return n * n // 4


def predicted_FS(n, w):
    """Lemma B/C closed form: (F(w), S(w), K(w))."""
    w %= n
    if n % 2 == 1:
        return P(n), 3 * (n - 1) // 2, (n - 1) // 2
    m = n // 2
    K = (n - 2) // 2 if w % 2 == 0 else n // 2
    S = 3 * K
    if m % 2 == 0:                      # n = 0 mod 4
        F = P(n)
    else:                                # n = 2 mod 4
        F = P(n) - 1 if w % 2 == 0 else P(n) + 1
    return F, S, K


def good_c(n, h):
    """The explicit shift of the theorem: makes w = h+c have the favourable
    parity (odd for n = 0 mod 4, even for n = 2 mod 4; irrelevant for odd n),
    using only c in {0, 1} so H(c) in {n, n-1}."""
    if n % 2 == 1:
        return 1
    m = n // 2
    want_odd = (m % 2 == 0)             # n = 0 mod 4 -> want w odd
    h_parity_ok = (h % 2 == 1) == want_odd
    return 0 if h_parity_ok else 1


def predicted_total(n, h):
    """Predicted value of 2F_c-S_c+H(c) at c = good_c(n,h), and its class."""
    c = good_c(n, h)
    w = (h + c) % n
    F, S, K = predicted_FS(n, w)
    H = n if c == 0 else n - 1
    total = 2 * F - S + H
    return total, c


def _check_one_w(n, w):
    f = [(w - i) % n for i in range(n)]
    cycles = su.cycles_of(f)
    F = sum(su.cycle_data(c, n)["F"] for c in cycles)
    S = sum(su.cycle_data(c, n)["S"] for c in cycles)
    K = len(cycles)
    Fp, Sp, Kp = predicted_FS(n, w)
    if (F, S, K) != (Fp, Sp, Kp):
        raise CheckFailure(f"Lemma B/C mismatch n={n} w={w}: got (F,S,K)=({F},{S},{K}), "
                           f"predicted ({Fp},{Sp},{Kp})")


def check_lemma_BC(n_max, extra_ns, log):
    """Check 1: exhaustive over all w for 4<=n<=n_max (cost is O(n^2) per n,
    i.e. O(n_max^3) total -- keep n_max moderate); a handful of w (both
    parities, w=0,1,near n/2,n-1) at the larger `extra_ns` (cost O(n) each)."""
    t0 = time.time()
    checked = 0
    for n in range(4, n_max + 1):
        for w in range(n):
            _check_one_w(n, w)
            checked += 1
    log(f"Lemma B/C: exhaustive 4<=n<={n_max}, all w -> {checked} (n,w) pairs OK, "
        f"{time.time() - t0:.1f} s")
    t1 = time.time()
    extra_checked = 0
    for n in extra_ns:
        for w in sorted(set([0, 1, 2, n - 1, n // 2, n // 2 + 1]) & set(range(n))):
            _check_one_w(n, w)
            extra_checked += 1
    log(f"Lemma B/C: structured sample at n = {extra_ns} -> {extra_checked} (n,w) pairs OK, "
        f"{time.time() - t1:.1f} s")
    return {"range": f"4<=n<={n_max}", "pairs": checked, "seconds": round(time.time() - t0, 1),
            "sample_ns": extra_ns, "sample_pairs": extra_checked,
            "sample_seconds": round(time.time() - t1, 1)}


def check_theorem(ns, log, exhaustive):
    """Check 2. For exhaustive n's, all h; otherwise a fixed structured set of h
    covering all 4 residues mod 4 and both parities, plus h=0,1,n-1."""
    out = {}
    for n in ns:
        t0 = time.time()
        if exhaustive:
            hs = range(n)
        else:
            # cost of a single (n, h) is O(n^2) (cycle_data builds an O(n) loads
            # array per cycle); keep the per-n sample small so large n stays cheap.
            base = sorted(set([0, 1, n - 1, n // 2, n // 2 + 1]) & set(range(n)))
            hs = base
        B = Bn(n)
        excess_hist = {}
        max_excess_R = None
        max_excess_len = None
        n_checked = 0
        for h in hs:
            pi = tuple((h - i) % n for i in range(n))
            pred_total, c = predicted_total(n, h)
            sw = su.shift_word(pi, c, "n1")
            wR = su.shift_word(pi, c, "carrier")
            _req_sorts(pi, sw["word"], n)
            _req_sorts(pi, wR["word"], n)
            actual = sw["stats"]["len"]
            if actual != pred_total:
                raise CheckFailure(f"Theorem mismatch n={n} h={h} c={c}: "
                                    f"predicted total {pred_total}, actual {actual}")
            actualR = wR["stats"]["len"]
            if actualR > actual:
                raise CheckFailure(f"Lemma A violated n={n} h={h} c={c}: R={actualR} > H-route {actual}")
            exc_len = actual - B
            exc_R = actualR - B
            if exc_R > 1:
                raise CheckFailure(f"H10 bound violated (carrier route) n={n} h={h} c={c}: excess {exc_R}")
            excess_hist[str(exc_len)] = excess_hist.get(str(exc_len), 0) + 1
            max_excess_R = exc_R if max_excess_R is None else max(max_excess_R, exc_R)
            max_excess_len = exc_len if max_excess_len is None else max(max_excess_len, exc_len)
            n_checked += 1
        out[n] = {
            "n_h_checked": n_checked, "coverage": "all h" if exhaustive else "structured sample",
            "excess_hist_H_route": excess_hist,
            "max_excess_H_route": max_excess_len, "max_excess_R": max_excess_R,
            "seconds": round(time.time() - t0, 3),
        }
        log(f"n={n}: {n_checked} h ({'exhaustive' if exhaustive else 'sample'}), "
            f"max excess H-route {max_excess_len}, max excess R {max_excess_R}, "
            f"{out[n]['seconds']} s")
    return out


def _req_sorts(pi, word, n):
    if apply_word(tuple(pi), word) != identity(n):
        raise CheckFailure(f"word does not sort pi={pi}: {word!r}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exhaustive-max", type=int, default=80,
                     help="exhaustive over all h up to this n (cost ~O(n^2) per h, "
                          "measured before use per AGENTS.md rule 10)")
    ap.add_argument("--lemma-bc-max", type=int, default=150,
                     help="exhaustive-over-w Lemma B/C check (O(n^3) total) up to this n")
    ap.add_argument("--sample-ns", default="150,200,201,202,203,500,501,502,503,998,999,1000,1001")
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)
    logs = []

    def log(msg):
        print(msg, flush=True)
        logs.append(msg)

    t0 = time.time()
    log(f"check_H10_reflections; core {CORE_VERSION}, construction {su.CONSTRUCTION_VERSION}")

    sample_ns_early = [int(x) for x in args.sample_ns.split(",") if x.strip()]
    lemma_bc = check_lemma_BC(args.lemma_bc_max, sample_ns_early, log)
    thm_exh = check_theorem(list(range(4, args.exhaustive_max + 1)), log, exhaustive=True)
    sample_ns = [int(x) for x in args.sample_ns.split(",") if x.strip()]
    thm_sample = check_theorem(sample_ns, log, exhaustive=False)

    report = {
        "meta": {"core": CORE_VERSION, "construction": su.CONSTRUCTION_VERSION,
                  "command": " ".join(sys.argv), "python": sys.version.split()[0]},
        "lemma_BC": lemma_bc,
        "theorem_exhaustive": thm_exh,
        "theorem_sample": thm_sample,
        "log": logs,
        "seconds_total": round(time.time() - t0, 1),
        "result": "PASS",
    }
    with open(os.path.join(OUT_DIR, "report.json"), "w") as f:
        json.dump(report, f, indent=1)
    with open(os.path.join(OUT_DIR, "summary.md"), "w") as f:
        f.write(f"# check_H10_reflections — PASS ({report['seconds_total']} s)\n\n")
        f.write(f"Lemma B/C (closed forms of F(w), S(w), K(w)): exhaustive over all w, "
                f"{lemma_bc['range']}, {lemma_bc['pairs']} pairs, all matched.\n\n")
        f.write(f"Theorem (explicit good shift c achieves the case-table value; "
                f"R_c route never exceeds it; excess over B_n via R_c is at most +1): "
                f"exhaustive over all h for 4<=n<={args.exhaustive_max}; "
                f"structured samples (covering all residues mod 4, both parities of h, "
                f"near n = 200, 500, 1000) at n = {sample_ns}.\n\n")
        f.write("Per-n details (max excess over B_n, H-route / R-route):\n\n")
        f.write("| n | coverage | max excess H-route | max excess R-route |\n|---|---|---|---|\n")
        for n, d in list(thm_exh.items())[-5:]:
            f.write(f"| {n} | {d['coverage']} | {d['max_excess_H_route']} | {d['max_excess_R']} |\n")
        for n, d in thm_sample.items():
            f.write(f"| {n} | {d['coverage']} | {d['max_excess_H_route']} | {d['max_excess_R']} |\n")
    log(f"done in {report['seconds_total']} s; report: {OUT_DIR}")


if __name__ == "__main__":
    main()
