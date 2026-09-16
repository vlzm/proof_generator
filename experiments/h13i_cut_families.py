"""H13-I candidate proof routes: does a restricted family of cuts, or plain
averaging over all n^2 cuts, already certify I(pi) <= floor((n-1)^2/4)?

I(pi) = min over all n^2 double cuts (q, c) of inv(pi relabelled on the line
starting at position q+1, values shifted by c) -- see experiments/line_model.py
and experiments/line_profile.c for the definition and the exhaustive
confirmation that max_pi I(pi) = floor((n-1)^2/4) at 4 <= n <= 10 (C33).

This script checks three candidate routes to a PROOF of that inequality:

1. Uniform averaging over all n^2 cuts: if avg_{q,c} inv(pi_{q,c}) were always
   <= floor((n-1)^2/4), pigeonhole would finish the proof with no further
   argument.  docs/notes/h13_line_model.md already states this fails on
   sigma_n's average and on id; this script finds the *worst* average exactly
   (over all pi, exhaustively at n <= 8) and its argmax.

2. Two restricted 1-parameter families (n cuts instead of n^2): value-only
   rotation (q fixed, no position rotation) and the diagonal family q = c.
   Both are weaker than the full n^2 search and are checked for whether they
   already reach the floor((n-1)^2/4) bound on all pi -- a positive result
   would reduce H13-I search space enormously.

Both routes are checked exhaustively at 4 <= n <= 8 (rule 9 AGENTS.md: cheap
enough to run to 8 directly) and the worst permutations for each route are
reported so the failure is a concrete, minimized counterexample, not just a
percentage.

Usage: python3 experiments/h13i_cut_families.py --nmax 8
Output: data/runs/h13i_cut_families/report.json, report.md.
"""

import argparse
import itertools
import json
import os
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
from moves import CORE_VERSION  # noqa: E402

VERSION = "h13i_cut_families-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h13i_cut_families")


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def line(pi, q, c):
    n = len(pi)
    return tuple((pi[(q + 1 + j) % n] - (q + 1 - c)) % n for j in range(n))


def run(n, log):
    t0 = time.time()
    bound = (n - 1) ** 2 // 4
    worst_avg, worst_avg_pi = -1.0, None
    worst_I, worst_I_pi = -1, None
    worst_valonly, worst_valonly_pi = -1, None
    worst_diag, worst_diag_pi = -1, None
    n_bad_avg = n_bad_valonly = n_bad_diag = 0
    cnt = 0
    for pi in itertools.permutations(range(n)):
        cnt += 1
        tot = 0
        best_I = None
        best_valonly = None
        best_diag = None
        for q in range(n):
            for c in range(n):
                v = inversions(line(pi, q, c))
                tot += v
                if best_I is None or v < best_I:
                    best_I = v
                if q == n - 1 and (best_valonly is None or v < best_valonly):
                    best_valonly = v
                if q == c and (best_diag is None or v < best_diag):
                    best_diag = v
        avg = tot / (n * n)
        if avg > worst_avg:
            worst_avg, worst_avg_pi = avg, pi
        if best_I > worst_I:
            worst_I, worst_I_pi = best_I, pi
        if best_valonly > worst_valonly:
            worst_valonly, worst_valonly_pi = best_valonly, pi
        if best_diag > worst_diag:
            worst_diag, worst_diag_pi = best_diag, pi
        if avg > bound:
            n_bad_avg += 1
        if best_valonly > bound:
            n_bad_valonly += 1
        if best_diag > bound:
            n_bad_diag += 1
    row = {
        "n": n, "bound_floor_(n-1)^2/4": bound, "count": cnt,
        "max_I": worst_I, "argmax_I": list(worst_I_pi),
        "max_avg": worst_avg, "argmax_avg": list(worst_avg_pi),
        "n_pi_with_avg_over_bound": n_bad_avg,
        "max_value_only_min": worst_valonly, "argmax_value_only_min": list(worst_valonly_pi),
        "n_pi_value_only_over_bound": n_bad_valonly,
        "max_diagonal_min": worst_diag, "argmax_diagonal_min": list(worst_diag_pi),
        "n_pi_diagonal_over_bound": n_bad_diag,
        "seconds": round(time.time() - t0, 1),
    }
    log(f"n={n}: bound={bound}; max I={worst_I} at {worst_I_pi} (must equal bound: {worst_I == bound}); "
        f"max avg={worst_avg:.3f} at {worst_avg_pi} ({n_bad_avg}/{cnt} pi exceed bound on average -- "
        f"averaging route FAILS); value-only family: max min={worst_valonly} at {worst_valonly_pi} "
        f"({n_bad_valonly}/{cnt} exceed bound -- family FAILS iff >0); diagonal family (q=c): "
        f"max min={worst_diag} at {worst_diag_pi} ({n_bad_diag}/{cnt} exceed bound); {row['seconds']} s")
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=4)
    ap.add_argument("--nmax", type=int, default=8)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    logf = open(os.path.join(OUT, "report.log"), "a")

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + "\n")
        logf.flush()

    log(f"== {VERSION} {CORE_VERSION} args={vars(args)}")
    rows = []
    for n in range(args.nmin, args.nmax + 1):
        rows.append(run(n, log))
    with open(os.path.join(OUT, f"report_n{args.nmin}_{args.nmax}.json"), "w") as f:
        json.dump({"version": VERSION, "core": CORE_VERSION, "rows": rows}, f, indent=1)


if __name__ == "__main__":
    main()
