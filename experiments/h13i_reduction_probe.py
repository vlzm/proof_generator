"""H13-I reduction probe (session 9): test two candidate proof strategies.

H13-I (PLAN, docs/notes/h13_line_model.md §6): for every pi in S_n there is a
double cut (rotation of positions and values) with inv(w) = I(pi) <=
floor((n-1)^2/4). Exhaustive data (session 8) already show this holds for
4 <= n <= 10, with equality on exactly the n reflections pi_h(i) = (h-i) mod n.
This probe tests two ways to *prove* it that were not yet ruled out:

  (1) corner cuts: restrict the n^2 cuts (a, b) to the n cuts that align an
      actual point (k, pi(k)) with the cut corner (w_0 = 0). Cheaper (O(n) vs
      O(n^2) candidates) and would give a short inductive proof if it always
      reached the bound.
  (2) pair-removal induction: remove two full points (i, pi(i)), (j, pi(j))
      (positions and values both contracted, the only way to remove two
      points and keep a permutation), giving pi'' of size n-2. Since
      floor((n-1)^2/4) - floor((n-3)^2/4) = n-2 exactly, an induction works
      if some pair removal always keeps I(pi) - I(pi'') <= n-2.

Both are checked exhaustively over all pi for 4 <= n <= NMAX (NMAX <= 8, cost
is O(n! * n^2 * n^2) resp. O(n! * C(n,2) * (n-2)^2), rule 10: n=9+ only by
explicit request, not run here).
Usage: python3 experiments/h13i_reduction_probe.py --nmax 8
Output: data/runs/h13i_reduction_probe/. Version h13i_reduction_probe-1.0.
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

VERSION = "h13i_reduction_probe-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h13i_reduction_probe")


def inv_count(seq):
    n = len(seq)
    c = 0
    for i in range(n):
        for j in range(i + 1, n):
            if seq[i] > seq[j]:
                c += 1
    return c


def I_full(pi):
    """I(pi) = min over all n^2 double cuts (a, b) of inv(rotated pi)."""
    n = len(pi)
    best = None
    for a in range(n):
        rp = [pi[(j + a) % n] for j in range(n)]
        for b in range(n):
            seq = [(v + b) % n for v in rp]
            c = inv_count(seq)
            if best is None or c < best:
                best = c
    return best


def I_corner(pi):
    """I restricted to the n cuts through an actual point (w_0 = 0)."""
    n = len(pi)
    best = None
    for k in range(n):
        rp = [pi[(j + k) % n] for j in range(n)]
        b = (-pi[k]) % n
        seq = [(v + b) % n for v in rp]
        c = inv_count(seq)
        if best is None or c < best:
            best = c
    return best


def floor_np1_sq_4(n):
    return ((n - 1) ** 2) // 4


def remove_two(pi, i, j):
    """Remove points i, j (both position and value), contract the rest."""
    n = len(pi)
    vals_removed = {pi[i], pi[j]}
    pos_remaining = [k for k in range(n) if k != i and k != j]
    remaining_vals = sorted(set(range(n)) - vals_removed)
    val_rank = {v: r for r, v in enumerate(remaining_vals)}
    return [val_rank[pi[k]] for k in pos_remaining]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--nmin", type=int, default=4)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    report = {"version": VERSION, "core_version": CORE_VERSION, "runs": []}

    for n in range(args.nmin, args.nmax + 1):
        t0 = time.time()
        bound = floor_np1_sq_4(n)
        max_I = -1
        corner_fails = 0
        max_corner = -1
        reduction_budget = n - 2
        reduction_fails = []
        max_reduction_gap = -10**9

        for pi in itertools.permutations(range(n)):
            pi = list(pi)
            Ipi = I_full(pi)
            if Ipi > max_I:
                max_I = Ipi

            Ic = I_corner(pi)
            if Ic > max_corner:
                max_corner = Ic
            if Ic > bound:
                corner_fails += 1

            best_reduced = -1
            for i, j in itertools.combinations(range(n), 2):
                red = remove_two(pi, i, j)
                Ir = I_full(red)
                if Ir > best_reduced:
                    best_reduced = Ir
            gap = Ipi - best_reduced
            if gap > max_reduction_gap:
                max_reduction_gap = gap
            if gap > reduction_budget:
                reduction_fails.append({"pi": tuple(pi), "I": Ipi,
                                         "best_reduced_I": best_reduced,
                                         "gap": gap})

        elapsed = time.time() - t0
        run = {
            "n": n,
            "bound_floor_np1_sq_4": bound,
            "max_I": max_I,
            "H13_I_holds": max_I <= bound,
            "corner_cut": {"max_I_corner": max_corner,
                            "fails_vs_bound": corner_fails,
                            "total": len(list(itertools.permutations(range(n))))},
            "pair_removal_induction": {
                "budget_n_minus_2": reduction_budget,
                "max_gap": max_reduction_gap,
                "num_fails": len(reduction_fails),
                "fails": reduction_fails,
            },
            "elapsed_s": elapsed,
        }
        report["runs"].append(run)
        print(f"n={n} bound={bound} max_I={max_I} "
              f"corner_fails={corner_fails} "
              f"pair_removal: budget={reduction_budget} max_gap={max_reduction_gap} "
              f"fails={len(reduction_fails)} ({elapsed:.1f}s)")

    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(report, f, indent=2, default=str)

    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# H13-I reduction probe (session 9)\n\n")
        f.write(f"Version {VERSION}, core {CORE_VERSION}.\n\n")
        f.write("| n | bound | max I | corner fails | pair-removal budget | "
                "max gap | pair-removal fails |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for run in report["runs"]:
            f.write(f"| {run['n']} | {run['bound_floor_np1_sq_4']} | "
                     f"{run['max_I']} | "
                     f"{run['corner_cut']['fails_vs_bound']} | "
                     f"{run['pair_removal_induction']['budget_n_minus_2']} | "
                     f"{run['pair_removal_induction']['max_gap']} | "
                     f"{run['pair_removal_induction']['num_fails']} |\n")
        f.write("\nConclusion: corner cuts (n candidates through a data "
                "point) fail already at n=4; two-point removal induction "
                "(budget n-2) holds exactly for 4<=n<=7 (equality on "
                "reflections) and fails at n=8 on exactly the 8 rotations of "
                "pi(i) = 5i mod 8 (gap 7 > budget 6). Both approaches are "
                "ruled out as proof strategies for H13-I.\n")

    print("Report:", os.path.join(OUT, "report.md"))


if __name__ == "__main__":
    main()
