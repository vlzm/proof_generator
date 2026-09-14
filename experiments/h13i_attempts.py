"""H13-I proof attempts (session 9): four candidate simplifications of the
double-cut search, tested exhaustively against the target bound
floor((n-1)^2/4).  All four fail (counterexamples below n = 9), i.e. none of
them alone can replace the full min over all n^2 cuts (q, c) in a proof of
H13-I.  See docs/notes/h13_line_model.md §7 for the write-up.

Recall the definitions (line_model.py, line_profile.c): for a cut edge
{q, q+1} and shift c, the relabelled line is
    w_j = (pi[(q+1+j) mod n] - (q+1-c)) mod n,  j = 0..n-1
I(pi) = min over (q, c) of inv(w).  H13-I conjectures max_pi I(pi) =
floor((n-1)^2/4) (VERIFIED 4 <= n <= 10, C33/C35).

Candidates (a = q+1, b = q+1-c, so w_j = (pi[(a+j) mod n] - b) mod n):

1. fix_a0: only vary b (a = 0 fixed).  Tests whether a single position-cut
   with optimal value-shift already suffices.
2. cut_at_data: force the cut to pass exactly through a data point, i.e.
   b = pi[a] for the chosen a (so w_0 = 0); only n choices instead of n^2.
3. footrule: replace inv(w) by the Diaconis-Graham footrule surrogate
   F(a,b) = sum_i |((i-a) mod n) - ((pi[i]-b) mod n)|.  Since inv <= F
   pointwise (Diaconis-Graham 1977), min_{a,b} F(a,b) <= floor((n-1)^2/4)
   would prove H13-I; tests whether the (analytically more tractable)
   footrule minimum obeys the same bound.
4. clean_interval: does pi always admit a circular position-interval of
   size m (for some a, m) whose value set is itself a circular interval
   (i.e. a "clean" cut with zero cross-inversions, giving upper bound
   C(m,2) + C(n-m,2) with zero cross term)?  Only useful if a *balanced*
   such m exists; measures the best bound this strategy can give (or None
   if pi has no clean interval at all beyond the trivial m = 1).

Usage: python3 experiments/h13i_attempts.py --nmax 8
Version h13i_attempts-1.0.
"""

import argparse
import itertools
import json
import os
import time
from math import comb

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "data", "runs", "h13i_attempts")
VERSION = "h13i_attempts-1.0"


def inv(seq):
    n = len(seq)
    return sum(1 for i in range(n) for j in range(i + 1, n) if seq[i] > seq[j])


def fix_a0(pi, n):
    best = None
    for b in range(n):
        iv = inv([(v - b) % n for v in pi])
        if best is None or iv < best:
            best = iv
    return best


def cut_at_data(pi, n):
    best = None
    for a in range(n):
        b = pi[a]
        seq = [(pi[(a + j) % n] - b) % n for j in range(n)]
        iv = inv(seq)
        if best is None or iv < best:
            best = iv
    return best


def footrule(pi, n):
    best = None
    for a in range(n):
        u = [(i - a) % n for i in range(n)]
        for b in range(n):
            v = [(pi[i] - b) % n for i in range(n)]
            F = sum(abs(u[i] - v[i]) for i in range(n))
            if best is None or F < best:
                best = F
    return best


def clean_interval(pi, n):
    best = None
    for a in range(n):
        for m in range(1, n):
            valset = set(pi[(a + t) % n] for t in range(m))
            if any(valset == set((s + t) % n for t in range(m)) for s in range(n)):
                bound_here = comb(m, 2) + comb(n - m, 2)
                if best is None or bound_here < best:
                    best = bound_here
    return best


CANDIDATES = {
    "fix_a0": (fix_a0, 9),
    "cut_at_data": (cut_at_data, 9),
    "footrule": (footrule, 8),
    "clean_interval": (clean_interval, 8),
}


def run(n, log):
    bound = (n - 1) ** 2 // 4
    total = 1
    for k in range(2, n + 1):
        total *= k
    stats = {}
    for name, (fn, nmax) in CANDIDATES.items():
        if n > nmax:
            continue
        t0 = time.time()
        worst, worst_pi, exceed, none_ct = -1, None, 0, 0
        for pi in itertools.permutations(range(n)):
            v = fn(pi, n)
            if v is None:
                none_ct += 1
                continue
            if v > worst:
                worst, worst_pi = v, pi
            if v > bound:
                exceed += 1
        stats[name] = {"worst": worst, "argmax": list(worst_pi) if worst_pi else None,
                        "exceed_count": exceed, "none_count": none_ct, "total": total,
                        "seconds": round(time.time() - t0, 1)}
        log(f"n={n} {name}: bound={bound} worst={worst} at {worst_pi} "
            f"exceed={exceed}/{total} none={none_ct} ({time.time()-t0:.1f}s)")
    return {"n": n, "bound": bound, "total": total, "candidates": stats}


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

    log(f"=== {VERSION} nmin={args.nmin} nmax={args.nmax} ===")
    rows = [run(n, log) for n in range(args.nmin, args.nmax + 1)]
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump({"version": VERSION, "rows": rows}, f, indent=1)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write(f"# H13-I proof attempts ({VERSION})\n\n")
        f.write("Four candidate simplifications, each tested exhaustively against "
                "`floor((n-1)^2/4)`. All four fail (counterexamples appear by n = 9 "
                "at the latest); see docs/notes/h13_line_model.md §7.\n\n")
        f.write("| n | bound | fix_a0 worst | cut_at_data worst | footrule worst | clean_interval worst |\n")
        f.write("|---|---|---|---|---|---|\n")
        for row in rows:
            c = row["candidates"]
            def cell(name):
                return str(c[name]["worst"]) if name in c else "-"
            f.write(f"| {row['n']} | {row['bound']} | {cell('fix_a0')} | {cell('cut_at_data')} | "
                    f"{cell('footrule')} | {cell('clean_interval')} |\n")
    log("done")


if __name__ == "__main__":
    main()
