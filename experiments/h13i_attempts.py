"""H13-I proof attempts: exhaustive tests of three candidate reductions.

H13-I conjectures I(pi) = min_{q,c} inv(line(pi,q,c)) <= floor((n-1)^2/4) for
every pi in S_n, n >= 4 (docs/notes/h13_line_model.md §6).  This script tests
three candidate simplifications that would make H13-I tractable, all of which
turn out to fail (exhaustively, on the ranges below):

1. Single freedom in c only, with q fixed at an arbitrary constant (n-1, i.e.
   pi read in its natural order): does max_pi min_c inv(pi shifted by c) stay
   <= floor((n-1)^2/4)?  FAILS from n = 7 (q freedom is load-bearing, not just
   a convenience -- this sharpens the existing "averaging over c alone does
   not work" note into an explicit worst-case counterexample per n).
2. The "w_0 = 0" coupling: for each q, take the unique c that makes the first
   line element 0 (no per-q search).  FAILS already at n = 4.
3. Coordinate-ascent local search on F(q, c) = n(n-1)/2 - inv(line(pi,q,c))
   from a fixed start (0, 0): does every local maximum reached already satisfy
   the bound?  FAILS with growing frequency from n = 4 (F is not unimodal on
   the (q, c) torus), so a "local max is global" argument is not available
   for free and any local-search-based proof needs a smarter starting point
   or an explicit escape from bad local optima.

Exhaustive over all pi for 4 <= n <= NMAX (default 8, budgeted per AGENTS.md
rule 10: n = 8 is 40320 permutations x <= n^2 evaluations, a few seconds).

Usage: python3 experiments/h13i_attempts.py --nmax 8
Output: data/runs/h13i_attempts/report.md, report.json.  Version h13i_attempts-1.0.
"""

import argparse
import itertools
import json
import math
import os
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT_DIR = os.path.join(ROOT, "data", "runs", "h13i_attempts")


def inv_count(w):
    n = len(w)
    c = 0
    for i in range(n):
        wi = w[i]
        for j in range(i + 1, n):
            if wi > w[j]:
                c += 1
    return c


def line(pi, q, c, n):
    off = (q + 1 - c) % n
    return [(pi[(q + 1 + j) % n] - off) % n for j in range(n)]


def target(n):
    return ((n - 1) ** 2) // 4


def attempt_fixed_q(pi, n, q_fixed):
    best = None
    for c in range(n):
        iv = inv_count(line(pi, q_fixed, c, n))
        if best is None or iv < best:
            best = iv
    return best


def attempt_w0_coupling(pi, n):
    best = None
    for q in range(n):
        c = (pi[(q + 1) % n] - 1) % n
        iv = inv_count(line(pi, q, c, n))
        if best is None or iv < best:
            best = iv
    return best


def local_search(pi, n, q0, c0):
    def F(q, c):
        return n * (n - 1) // 2 - inv_count(line(pi, q, c, n))

    q, c = q0, c0
    cur = F(q, c)
    improved = True
    while improved:
        improved = False
        for nq, nc in ((q + 1) % n, c), ((q - 1) % n, c), (q, (c + 1) % n), (q, (c - 1) % n):
            v = F(nq, nc)
            if v > cur:
                cur, q, c = v, nq, nc
                improved = True
                break
    return n * (n - 1) // 2 - cur


def run(nmax):
    rows = []
    for n in range(4, nmax + 1):
        t = target(n)
        total = math.factorial(n)
        r1_fails = r1_worst = 0
        r2_fails = r2_worst = 0
        r3_fails = r3_worst = 0
        r1_example = r2_example = r3_example = None
        t0 = time.time()
        for pi in itertools.permutations(range(n)):
            v1 = attempt_fixed_q(pi, n, n - 1)
            if v1 > t:
                r1_fails += 1
                if v1 - t > r1_worst:
                    r1_worst = v1 - t
                    r1_example = pi
            v2 = attempt_w0_coupling(pi, n)
            if v2 > t:
                r2_fails += 1
                if v2 - t > r2_worst:
                    r2_worst = v2 - t
                    r2_example = pi
            v3 = local_search(pi, n, 0, 0)
            if v3 > t:
                r3_fails += 1
                if v3 - t > r3_worst:
                    r3_worst = v3 - t
                    r3_example = pi
        rows.append(dict(
            n=n, target=t, total=total, elapsed=time.time() - t0,
            fixed_q=dict(fails=r1_fails, worst_gap=r1_worst, example=r1_example),
            w0_coupling=dict(fails=r2_fails, worst_gap=r2_worst, example=r2_example),
            local_search=dict(fails=r3_fails, worst_gap=r3_worst, example=r3_example),
        ))
        print(n, rows[-1])
    return rows


def write_report(rows, nmax):
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(os.path.join(OUT_DIR, "report.json"), "w") as f:
        json.dump(rows, f, indent=2)
    lines = []
    lines.append("# H13-I proof attempts: report\n")
    lines.append(f"h13i_attempts-1.0. Exhaustive 4 <= n <= {nmax}. All pi.\n")
    lines.append("| n | target | fixed_q fails | fixed_q worst gap | example | "
                  "w0 fails | w0 worst gap | example | local-search fails | "
                  "local-search worst gap | example |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        lines.append(
            f"| {r['n']} | {r['target']} | {r['fixed_q']['fails']} | "
            f"{r['fixed_q']['worst_gap']} | {r['fixed_q']['example']} | "
            f"{r['w0_coupling']['fails']} | {r['w0_coupling']['worst_gap']} | "
            f"{r['w0_coupling']['example']} | {r['local_search']['fails']} | "
            f"{r['local_search']['worst_gap']} | {r['local_search']['example']} |"
        )
    lines.append("")
    lines.append("Conclusion: all three candidate reductions of H13-I fail exhaustively "
                  "from small n (see docs/notes/h13_line_model.md §7 for the verdict). "
                  "The full 2-parameter search over (q, c) remains necessary and no "
                  "closed-form / single-pass construction was found this session.")
    with open(os.path.join(OUT_DIR, "report.md"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8)
    args = ap.parse_args()
    rows = run(args.nmax)
    write_report(rows, args.nmax)
