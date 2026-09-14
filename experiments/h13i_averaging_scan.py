"""H13-I, session 9: how badly do averaging-based candidate proofs of H13-I
fail, and how far from the target do they fail?  (docs/notes/h13i_session9.md,
docs/notes/h13_line_model.md Section 6, docs/proofs/C37_quadrant_pair_identity.md)

Two averaging candidates, both ruled out here quantitatively (not just "fails
for a special family" as already known from h13_line_model.md):

  (A) full 2D average: avg_{q,c} inv(w_{q,c}) <= floor((n-1)^2/4)?
      Computed exactly via the C37 pair identity (sum over C(n,2) pairs of
      a(n-b)+(n-a)b, divided by n^2) -- no need to build the n^2 lines.
  (B) best-q, then average over c: min_q [ (1/n) sum_c inv(w_{q,c}) ]
      <= floor((n-1)^2/4)?  Built directly (n^2 lines per pi; only feasible
      exhaustively for small n, so n = 9, 10 are a random sample here with
      the sample size logged -- rule 9).

Both are much weaker requirements than H13-I itself (they ask for an
averaging argument to already prove the bound; H13-I only needs existence of
one good cut, already known VERIFIED 4 <= n <= 10 via exhaustive BFS/search
in data/runs/line_profile).  Both fail on almost every permutation already at
small n, not just on reflections or antipodal transpositions.

Note (found while running this script, proved in docs/proofs/
C37_quadrant_pair_identity.md Lemma 2): (A) and (B) are not just numerically
close here -- they are *exactly* equal for every pi, because sum_c
inv(w_{q,c}) does not depend on q at all (row sums of the (q,c) grid are all
equal).  So (B) has no real extra freedom over (A); the table below shows
this as identical fail-counts and worst-case values, not a coincidence of
the sample.

Usage: python3 experiments/h13i_averaging_scan.py [--nmax 9] [--sample-n 9,10] [--samples 20000] [--seed 1]
Output: data/runs/h13i_averaging_scan/report.json, report.md.
"""

import argparse
import itertools
import json
import os
import random
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
VERSION = "h13i_averaging_scan-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h13i_averaging_scan")


def bound(n):
    return ((n - 1) ** 2) // 4


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def full_2d_average(pi):
    """avg_{q,c} inv(w_{q,c}), via the C37 pair identity (exact, O(n^2))."""
    n = len(pi)
    total = 0
    for i1 in range(n):
        for i2 in range(i1 + 1, n):
            a = (i2 - i1) % n
            b = (pi[i2] - pi[i1]) % n
            total += a * (n - b) + (n - a) * b
    return total / (n * n)


def min_q_avg_c(pi):
    """min over q of (1/n) sum_c inv(w_{q,c})."""
    n = len(pi)
    best = None
    for q in range(n):
        u = [pi[(q + 1 + j) % n] for j in range(n)]
        total = 0
        for c in range(n):
            shift = (q + 1 - c) % n
            w = [(v - shift) % n for v in u]
            total += inversions(w)
        avg = total / n
        if best is None or avg < best:
            best = avg
    return best


def exhaustive_part(nmax, log):
    rows = []
    for n in range(4, nmax + 1):
        b = bound(n)
        t0 = time.time()
        total = 0
        fail_a = 0
        fail_b = 0
        worst_a = -1.0
        worst_b = -1.0
        worst_a_pi = None
        worst_b_pi = None
        for pi in itertools.permutations(range(n)):
            total += 1
            va = full_2d_average(pi)
            if va > b + 1e-9:
                fail_a += 1
            if va > worst_a:
                worst_a, worst_a_pi = va, pi
            vb = min_q_avg_c(pi)
            if vb > b + 1e-9:
                fail_b += 1
            if vb > worst_b:
                worst_b, worst_b_pi = vb, pi
        dt = time.time() - t0
        row = {
            "n": n, "bound": b, "total": total,
            "fail_full2d": fail_a, "worst_full2d": worst_a, "worst_full2d_pi": worst_a_pi,
            "fail_bestq_avgc": fail_b, "worst_bestq_avgc": worst_b, "worst_bestq_avgc_pi": worst_b_pi,
            "seconds": dt, "coverage": "exhaustive",
        }
        rows.append(row)
        log(f"n={n} bound={b}: (A) full-2D-avg fails on {fail_a}/{total} "
            f"(worst {worst_a:.3f} at {worst_a_pi}); "
            f"(B) best-q avg-c fails on {fail_b}/{total} "
            f"(worst {worst_b:.3f} at {worst_b_pi})  [{dt:.1f}s]")
    return rows


def sample_part(ns, samples, seed, log):
    rows = []
    rng = random.Random(seed)
    for n in ns:
        b = bound(n)
        t0 = time.time()
        fail_a = 0
        fail_b = 0
        worst_a = -1.0
        worst_b = -1.0
        for _ in range(samples):
            pi = list(range(n))
            rng.shuffle(pi)
            pi = tuple(pi)
            va = full_2d_average(pi)
            if va > b + 1e-9:
                fail_a += 1
            worst_a = max(worst_a, va)
            vb = min_q_avg_c(pi)
            if vb > b + 1e-9:
                fail_b += 1
            worst_b = max(worst_b, vb)
        dt = time.time() - t0
        row = {
            "n": n, "bound": b, "total": samples,
            "fail_full2d": fail_a, "worst_full2d": worst_a,
            "fail_bestq_avgc": fail_b, "worst_bestq_avgc": worst_b,
            "seconds": dt, "coverage": f"random sample, {samples} uniform pi (seed={seed})",
        }
        rows.append(row)
        log(f"n={n} bound={b} [SAMPLE {samples}]: (A) fails on {fail_a}/{samples} "
            f"(worst {worst_a:.3f}); (B) fails on {fail_b}/{samples} "
            f"(worst {worst_b:.3f})  [{dt:.1f}s]")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=9)
    ap.add_argument("--sample-n", type=str, default="10")
    ap.add_argument("--samples", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    lines = []

    def log(s):
        print(s)
        lines.append(s)

    rows = exhaustive_part(args.nmax, log)
    sample_ns = [int(x) for x in args.sample_n.split(",") if x.strip()]
    sample_rows = sample_part(sample_ns, args.samples, args.seed, log)

    report = {
        "version": VERSION,
        "nmax": args.nmax,
        "sample_n": sample_ns,
        "samples": args.samples,
        "seed": args.seed,
        "exhaustive": rows,
        "sampled": sample_rows,
        "log": lines,
    }
    with open(os.path.join(OUT, "report.json"), "w") as fh:
        json.dump(report, fh, indent=2)

    with open(os.path.join(OUT, "report.md"), "w") as fh:
        fh.write(f"# h13i_averaging_scan report ({VERSION})\n\n")
        fh.write("Two averaging candidates for H13-I, both ruled out quantitatively.\n\n")
        fh.write("## Exhaustive (4 <= n <= {})\n\n".format(args.nmax))
        fh.write("| n | bound | total | (A) full-2D-avg fails | worst (A) | "
                  "(B) best-q avg-c fails | worst (B) |\n")
        fh.write("|---|---|---|---|---|---|---|\n")
        for r in rows:
            fh.write(f"| {r['n']} | {r['bound']} | {r['total']} | "
                      f"{r['fail_full2d']} | {r['worst_full2d']:.3f} | "
                      f"{r['fail_bestq_avgc']} | {r['worst_bestq_avgc']:.3f} |\n")
        fh.write("\n## Sampled\n\n")
        fh.write("| n | bound | samples | (A) fails | worst (A) | (B) fails | worst (B) |\n")
        fh.write("|---|---|---|---|---|---|---|\n")
        for r in sample_rows:
            fh.write(f"| {r['n']} | {r['bound']} | {r['total']} | "
                      f"{r['fail_full2d']} | {r['worst_full2d']:.3f} | "
                      f"{r['fail_bestq_avgc']} | {r['worst_bestq_avgc']:.3f} |\n")

    print(f"\nReport: {os.path.join(OUT, 'report.md')}")


if __name__ == "__main__":
    main()
