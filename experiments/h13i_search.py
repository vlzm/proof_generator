"""H13-I: does the gap between I(pi) and the "fix position-cut a = 0, vary
only the value-cut b" restriction stay bounded as n grows?

Companion to experiments/h13i_coupling.py, which showed the gap is already
+1 at n = 7, 8, 9 (exhaustive).  This script checks larger n two ways:

1. Exhaustive over reflections pi_h(i) = (h - i) mod n (the known I(pi)
   extremizers, C33/C35): is the fixed-a0 restriction exact on them?
2. Randomized hill-climbing (2-opt on the permutation, fixed seed) trying to
   maximize (fixed_a0_min(pi) - floor((n-1)^2/4)): does the reachable gap
   grow with n, or plateau?

This is a search for a growing counterexample to "one free parameter (b,
with position cut fixed at 0) suffices for H13-I", not a proof of anything;
hill-climbing only demonstrates a lower bound on the true worst-case gap for
this restricted family.

Usage: python3 experiments/h13i_search.py --iters 4000 --seed 2
Output: appended to data/runs/h13i_coupling/report.md.
Version h13i_search-1.0 (session 9).
"""

import argparse
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h13i_coupling import fixed_a0_min  # noqa: E402

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def reflection(n, h):
    return [(h - i) % n for i in range(n)]


def reflections_scan(ns):
    rows = []
    for n in ns:
        bound = (n - 1) ** 2 // 4
        worst_gap = None
        worst_h = None
        for h in range(n):
            p = reflection(n, h)
            gap = fixed_a0_min(p, n) - bound
            if worst_gap is None or gap > worst_gap:
                worst_gap, worst_h = gap, h
        rows.append((n, bound, worst_gap, worst_h))
    return rows


def hill_climb(n, iters, rng):
    bound = (n - 1) ** 2 // 4
    curp = list(range(n))[::-1]
    cur = fixed_a0_min(curp, n)
    best, bestp = cur, curp[:]
    for _ in range(iters):
        i, j = rng.sample(range(n), 2)
        curp[i], curp[j] = curp[j], curp[i]
        v = fixed_a0_min(curp, n)
        if v >= cur:
            cur = v
        else:
            curp[i], curp[j] = curp[j], curp[i]
        if cur > best:
            best, bestp = cur, curp[:]
    return best - bound, bestp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iters", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=2)
    ap.add_argument("--reflection_ns", type=str,
                     default="10,15,20,30,40,50,60,80,100")
    ap.add_argument("--climb_ns", type=str, default="9,12,15,20,25,30")
    args = ap.parse_args()

    refl_ns = [int(x) for x in args.reflection_ns.split(",")]
    climb_ns = [int(x) for x in args.climb_ns.split(",")]

    t0 = time.time()
    refl_rows = reflections_scan(refl_ns)

    rng = random.Random(args.seed)
    climb_rows = []
    for n in climb_ns:
        gap, bestp = hill_climb(n, args.iters, rng)
        climb_rows.append((n, gap, tuple(bestp)))
    elapsed = time.time() - t0

    out_dir = os.path.join(ROOT, "data", "runs", "h13i_coupling")
    os.makedirs(out_dir, exist_ok=True)
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    lines = [f"\n## Search run {ts} ({elapsed:.1f} s), seed={args.seed}, iters={args.iters}\n"]
    lines.append("\n### Reflections: fixed-a0 gap vs floor((n-1)^2/4) (exhaustive over h)\n")
    for n, bound, gap, h in refl_rows:
        lines.append(f"- n={n}: bound={bound}, worst gap={gap} at h={h}\n")
    lines.append("\n### Hill-climb (2-opt from reversed identity, "
                  f"{args.iters} moves): best gap found\n")
    for n, gap, bestp in climb_rows:
        lines.append(f"- n={n}: gap={gap}, example={bestp}\n")

    report_path = os.path.join(out_dir, "report.md")
    with open(report_path, "a") as f:
        f.writelines(lines)
    print("".join(lines))
    print(f"appended to {report_path}")


if __name__ == "__main__":
    main()
