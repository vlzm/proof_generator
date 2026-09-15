"""Toric double-cut search (H13-I, session 9): a constructive candidate cut.

H13-I conjectures: for every pi in S_n there is a double cut (q, c) (PROBLEM /
docs/notes/h13_line_model.md S0: line w_j = pi(q+1+j) - (q+1-c) mod n) with
I(pi) = min_{q,c} inv(w) <= floor((n-1)^2/4).  Exhaustive search over all n^2
cuts already VERIFIED this at 4 <= n <= 10 (docs/notes/h13_line_model.md S1.1,
experiments/line_profile.c); no proof for general n is known.  Prior narrow
cut families were REFUTED (docs/notes/h13_line_model.md S1.7): induction on
one element, and the n-cut family "first line element has value t".

This module tests a new, richer candidate family, aiming at a constructive
proof rather than another exhaustive table.

1. Algebraic identity (elementary; proved in docs/notes/h13_line_model.md S7):
   for a fixed linear sequence u (a permutation of 0..n-1) and value-rotations
   w(t)_j = (u_j - t) mod n, t = 0..n-1,

       (1/n) * sum_t inv(w(t))  =  inv(u) + Q(u) / n,   where
       Q(u) = sum_j (2j - (n-1)) * u_j.

   So min_t inv(w(t)) <= inv(u) + Q(u) / n for every u (average >= min).

2. Candidate positions q for the OTHER cut: let e_i = 2*pi(i) - (n-1) (i =
   0..n-1, sum e_i = 0) and P_k = sum_{i<k} e_i (k = 0..n, P_0 = P_n = 0).  A
   cyclic *peak* is an index k in {0,...,n-1} with P_k >= P_{k-1} and
   P_k >= P_{k+1} (indices mod n); equivalently pi((k-1) mod n) >= (n-1)/2 >=
   pi(k) (the raw sequence crosses the median going down right at the cut).
   Candidate cuts: q = k - 1 for every peak k (there is always at least one:
   the global argmax of P is always a peak).

Claim tested (call it H13-I''): among the peak cuts, taking for each the best
value-shift via direct search over t (not just the averaging bound above),
min over peaks k of min_t inv(w) <= floor((n-1)^2/4).

Result of this session's exhaustive search: TRUE for every pi at 4 <= n <=
10 (all permutations, all peaks, all n choices of t: still nowhere near the
full n^2 search since peaks are typically few, but this is a *sufficiency*
check of a specific rule, not a proof).  NOT proved for general n: the global
argmax of P alone is insufficient (found counterexamples at n = 7, 9, e.g.
pi = (0,6,5,4,3,2,8,1,7) at n = 9 needs a non-global peak, P = -2, while both
global peaks P = 0 give inv = 17 > floor(8^2/4) = 16); no argument found this
session for why some peak always works when several disagree, and the
averaging identity above is not tight enough by itself (average over q AND t
jointly reduces to the known-insufficient full (q,c) average, docs/notes/
h13_line_model.md S6).

Usage: python3 experiments/toric_cut_peak.py --nmin 4 --nmax 10
Output: data/runs/toric_cut_peak/report.json, report.md.
Version toric_cut_peak-1.0.
"""

import argparse
import itertools
import json
import os
import time
from fractions import Fraction

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
VERSION = "toric_cut_peak-1.0"
OUT = os.path.join(ROOT, "data", "runs", "toric_cut_peak")


def inv_count(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def best_t_inv(u):
    """min_t inv((u - t) mod n), t = 0..n-1, via O(1) incremental update."""
    n = len(u)
    iv = inv_count(u)
    best = iv
    pos_of = [0] * n
    for idx, val in enumerate(u):
        pos_of[val] = idx
    for t in range(n - 1):
        j0 = pos_of[t]
        iv += n - 1 - 2 * j0
        if iv < best:
            best = iv
    return best


def avg_t_inv(u):
    n = len(u)
    Q = sum((2 * j - (n - 1)) * u[j] for j in range(n))
    return inv_count(u) + Fraction(Q, n)


def peaks_of(pi):
    n = len(pi)
    e = [2 * x - (n - 1) for x in pi]
    P = [0] * (n + 1)
    for k in range(n):
        P[k + 1] = P[k] + e[k]
    assert P[n] == 0
    out = []
    for k in range(n):
        pk = P[k]
        pprev = P[k - 1] if k > 0 else P[n - 1]
        pnext = P[k + 1] if k < n else P[0]
        if pk >= pprev and pk >= pnext:
            out.append(k)
    return out, P


def rotate(pi, k):
    n = len(pi)
    return [pi[(k + j) % n] for j in range(n)]


def peak_construction_best(pi):
    """min over peaks k of best_t_inv(rotate(pi, k)); also report per-peak values."""
    n = len(pi)
    peaks, P = peaks_of(pi)
    vals = [(k, best_t_inv(rotate(pi, k))) for k in peaks]
    best = min(v for _, v in vals)
    global_max = max(P[:n])
    global_argmax_vals = [v for k, v in vals if P[k] == global_max]
    return best, vals, global_argmax_vals


def full_search_I(pi):
    """Reference (slow, O(n^3)): true I(pi) over all n^2 cuts."""
    n = len(pi)
    best = None
    for q in range(n):
        for c in range(n):
            shift = (q + 1 - c) % n
            w = [(pi[(q + 1 + j) % n] - shift) % n for j in range(n)]
            iv = inv_count(w)
            if best is None or iv < best:
                best = iv
    return best


def run(nmin, nmax, cross_check_full, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    per_n = []
    counterexamples_global_argmax = []
    for n in range(nmin, nmax + 1):
        bound = ((n - 1) ** 2) // 4
        t0 = time.time()
        total = 0
        worst_peak = -1
        fails_peak = 0
        fails_global_argmax_only = 0
        worst_pi = None
        cross_check_ok = True
        for pi in itertools.permutations(range(n)):
            total += 1
            pi = list(pi)
            best, vals, global_vals = peak_construction_best(pi)
            if best > worst_peak:
                worst_peak = best
                worst_pi = tuple(pi)
            if best > bound:
                fails_peak += 1
            if min(global_vals) > bound:
                fails_global_argmax_only += 1
                if len(counterexamples_global_argmax) < 3:
                    counterexamples_global_argmax.append(
                        {"n": n, "pi": tuple(pi), "global_argmax_vals": global_vals,
                         "all_peak_vals": vals, "bound": bound}
                    )
            if cross_check_full and cross_check_ok:
                if full_search_I(pi) > bound:
                    cross_check_ok = False
        dt = time.time() - t0
        per_n.append({
            "n": n, "bound": bound, "total_permutations": total,
            "worst_peak_construction": worst_peak, "worst_pi": worst_pi,
            "fails_peak_construction": fails_peak,
            "fails_global_argmax_only": fails_global_argmax_only,
            "cross_check_full_I_le_bound": cross_check_ok if cross_check_full else None,
            "time_seconds": round(dt, 3),
        })
        print(per_n[-1])
    failing_ns = [str(r["n"]) for r in per_n if r["fails_global_argmax_only"] > 0]
    report = {
        "version": VERSION,
        "claim": "H13-I: I(pi) <= floor((n-1)^2/4) for all pi in S_n",
        "candidate_rule": "peaks of prefix sum of e_i=2pi(i)-(n-1); best value-shift per peak",
        "per_n": per_n,
        "counterexamples_global_argmax_only": counterexamples_global_argmax,
        "conclusion": (
            "peak-construction (all peaks, best of each) matches the bound "
            "exactly (no failures) for every pi at %d<=n<=%d; the "
            "restricted global-argmax-only version fails at n in {%s} "
            "(see counterexamples); no general-n proof found this session"
            % (nmin, nmax, ", ".join(failing_ns) if failing_ns else "none observed")
        ),
    }
    with open(os.path.join(out_dir, "report.json"), "w") as f:
        json.dump(report, f, indent=2, default=str)
    with open(os.path.join(out_dir, "report.md"), "w") as f:
        f.write("# toric_cut_peak: peak-rotation candidate for H13-I\n\n")
        f.write("Version: %s\n\n" % VERSION)
        f.write("| n | bound floor((n-1)^2/4) | worst (peak constr.) | fails | fails (global argmax only) | permutations | time (s) |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for r in per_n:
            f.write("| %d | %d | %d | %d | %d | %d | %.2f |\n" % (
                r["n"], r["bound"], r["worst_peak_construction"], r["fails_peak_construction"],
                r["fails_global_argmax_only"], r["total_permutations"], r["time_seconds"]))
        f.write("\n" + report["conclusion"] + "\n\n")
        if counterexamples_global_argmax:
            f.write("## Counterexamples to the global-argmax-only restriction\n\n")
            for ce in counterexamples_global_argmax:
                f.write("- n=%d, pi=%s, bound=%d: global-argmax peaks give %s, "
                        "but the full peak set gives %s\n" % (
                            ce["n"], ce["pi"], ce["bound"], ce["global_argmax_vals"],
                            [v for _, v in ce["all_peak_vals"]]))
    return report


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=4)
    ap.add_argument("--nmax", type=int, default=9)
    ap.add_argument("--cross-check-full", action="store_true",
                     help="also run the O(n^3)-per-pi full n^2 cut search as an "
                          "independent check of I(pi) <= bound (slow beyond n~8)")
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()
    run(args.nmin, args.nmax, args.cross_check_full, args.out)
