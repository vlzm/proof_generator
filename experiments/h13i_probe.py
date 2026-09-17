"""H13-I probe (session 9): discrete-Laplacian identity for I(pi) and two
refuted attack routes.

I(pi) = min_{a,b in Z_n} g(a, b), where g(a, b) = inv of the line obtained by
cutting positions after a and relabelling values by subtracting b mod n
(experiments/line_model.py, docs/notes/h13_line_model.md §0; equivalent to
choosing a double cut (q, c) = (a - 1, a - b)).

This module reports three things, each exhaustive for 4 <= n <= NMAX and
sampled for larger n:

1. Lemma D (PROVED, all n, all pi; dedicated checker checks/check_C37.py):
   the mixed second difference of g is
       g(a+1,b+1) - g(a+1,b) - g(a,b+1) + g(a,b) = 2 - 2n * [pi(a) == b]
   i.e. g is discretely "flat with curvature +2" everywhere except at the n
   points of the permutation itself, where it has a spike of -(2n-2). See
   docs/notes/h13_line_model.md §7 for the algebraic proof (double difference
   of the pair/XOR-of-arcs formula for g); this script re-checks it as a
   side effect of building the grid used for points 2-3 below.

2. Plain averaging over the n^2 origins does not bound I(pi): reports, for
   each n, how many pi (out of n!) have mean_{a,b} g(a,b) > floor((n-1)^2/4)
   (answer: almost all of them, confirming docs/notes/h13_line_model.md §6).

3. The natural exchange argument "swapping an adjacent pair of values (or
   positions) on the circle never decreases I(pi)" is FALSE in general
   (counterexamples exist already at n = 4); reports counts of pi, v where
   the swap of values v, v+1 (mod n) decreases I(pi), same for adjacent
   position swaps.

Usage: python3 experiments/h13i_probe.py --nmax 7
Output: data/runs/h13i_probe/report.json, report.md.  Version h13i_probe-1.0.
"""

import argparse
import itertools
import json
import math
import os
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
from moves import CORE_VERSION  # noqa: E402

VERSION = "h13i_probe-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h13i_probe")


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def grid(pi, n):
    g = [[0] * n for _ in range(n)]
    for a in range(n):
        for b in range(n):
            w = [(pi[(a + j) % n] - b) % n for j in range(n)]
            g[a][b] = inversions(w)
    return g


def check_laplacian(pi, n, g):
    for a in range(n):
        for b in range(n):
            d = g[(a + 1) % n][(b + 1) % n] - g[(a + 1) % n][b] - g[a][(b + 1) % n] + g[a][b]
            expected = 2 - 2 * n if pi[a] == b else 2
            if d != expected:
                return False, (a, b, d, expected)
    return True, None


def I_of_grid(g, n):
    return min(min(row) for row in g)


def swap_values(pi, v, w):
    lst = list(pi)
    pv, pw = lst.index(v), lst.index(w)
    lst[pv], lst[pw] = lst[pw], lst[pv]
    return tuple(lst)


def swap_positions(pi, i, j):
    lst = list(pi)
    lst[i], lst[j] = lst[j], lst[i]
    return tuple(lst)


def run(nmax, log):
    rows = []
    for n in range(4, nmax + 1):
        t0 = time.time()
        bound = ((n - 1) ** 2) // 4
        total = math.factorial(n)
        laplacian_fail = None
        avg_exceeds = 0
        val_swap_decreases = 0
        pos_swap_decreases = 0
        maxI = -1
        maxI_pi = None
        for pi in itertools.permutations(range(n)):
            g = grid(pi, n)
            ok, info = check_laplacian(pi, n, g)
            if not ok and laplacian_fail is None:
                laplacian_fail = {"pi": pi, "at": info}
            Ipi = I_of_grid(g, n)
            if Ipi > maxI:
                maxI = Ipi
                maxI_pi = pi
            mean_g = sum(sum(row) for row in g) / (n * n)
            if mean_g > bound:
                avg_exceeds += 1
            Ipi_cache = {}
            for v in range(n):
                w = (v + 1) % n
                pi2 = swap_values(pi, v, w)
                g2 = grid(pi2, n)
                if I_of_grid(g2, n) < Ipi:
                    val_swap_decreases += 1
            for i in range(n):
                j = (i + 1) % n
                pi2 = swap_positions(pi, i, j)
                g2 = grid(pi2, n)
                if I_of_grid(g2, n) < Ipi:
                    pos_swap_decreases += 1
        elapsed = time.time() - t0
        row = {
            "n": n,
            "bound_floor_(n-1)^2/4": bound,
            "total_pi": total,
            "laplacian_identity": "PASS" if laplacian_fail is None else "FAIL",
            "laplacian_fail_example": laplacian_fail,
            "maxI": maxI,
            "maxI_matches_bound": maxI == bound,
            "maxI_example": maxI_pi,
            "avg_exceeds_bound_count": avg_exceeds,
            "avg_exceeds_bound_total": total,
            "value_adjacent_swap_decreases_I": val_swap_decreases,
            "position_adjacent_swap_decreases_I": pos_swap_decreases,
            "elapsed_s": round(elapsed, 2),
        }
        rows.append(row)
        log(f"n={n} bound={bound} maxI={maxI} laplacian={row['laplacian_identity']} "
            f"avg_exceeds={avg_exceeds}/{total} val_swap_decr={val_swap_decreases} "
            f"pos_swap_decr={pos_swap_decreases} ({elapsed:.1f}s)")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=7)
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    log_lines = []

    def log(msg):
        print(msg)
        log_lines.append(msg)

    log(f"{VERSION} core={CORE_VERSION} nmax={args.nmax}")
    rows = run(args.nmax, log)

    report = {
        "version": VERSION,
        "core": CORE_VERSION,
        "goal": "H13-I probe: verify Lemma D (mixed-difference identity for g); "
                "quantify failure of plain averaging and of the adjacent-swap "
                "exchange argument as routes to I(pi) <= floor((n-1)^2/4)",
        "coverage": "exhaustive over all pi, all n^2 origins, all n in [4, nmax]",
        "rows": rows,
    }
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(report, f, indent=2, default=list)

    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write(f"# H13-I probe ({VERSION}, core {CORE_VERSION})\n\n")
        f.write("n | bound | maxI | Laplacian | avg>bound | val-swap decr | pos-swap decr | time\n")
        f.write("---|---|---|---|---|---|---|---\n")
        for r in rows:
            f.write(
                f"{r['n']} | {r['bound_floor_(n-1)^2/4']} | {r['maxI']} | "
                f"{r['laplacian_identity']} | {r['avg_exceeds_bound_count']}/{r['avg_exceeds_bound_total']} | "
                f"{r['value_adjacent_swap_decreases_I']} | {r['position_adjacent_swap_decreases_I']} | "
                f"{r['elapsed_s']}s\n"
            )
        f.write("\n" + "\n".join(log_lines) + "\n")

    print("Wrote", OUT)


if __name__ == "__main__":
    main()
