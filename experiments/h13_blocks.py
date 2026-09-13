"""H13 step (2), part A: does a min-N_X line scheme (N_X = I(pi), N_rot <= I(pi) + O(n))
exist for ALL pi, generalizing the C32 zigzag beyond reflections?

Counter-family: pi_k(n) = identity except k disjoint local transpositions
(swap positions {0,1}, {floor(n/k), floor(n/k)+1}, ..., {floor((k-1)n/k), +1}),
evenly spaced with one swap edge placed exactly at 0.  I(pi_k) = k for all
tested n (the swaps are pairwise non-interacting, each contributes exactly one
inversion on the best line).  Any word with N_X = I(pi_k) = k must place its
(unique, since fewer is impossible and the construction goal wants exactly the
minimum) swap at each of the k required edges; between two edges the head needs
at least the circular distance, and the start is fixed at head 0.  This is the
"shortest path from a fixed start visiting k points on a circle" problem: the
minimum is n - maxgap, where maxgap is the largest gap among the k points
(0 is one of them by construction, so no separate correction term).  For fixed
k >= 3 and n -> infinity, maxgap = n/k (evenly spaced), so the forced rotation
cost is n(1 - 1/k) + O(1), which grows strictly faster than floor(n/2) + O(1)
whenever k >= 3 (since 1 - 1/k > 1/2 iff k > 2).  This refutes the existence of
a bounded-overhead (N_rot <= I(pi) + floor(n/2) + O(1)) minimal-N_X scheme for
general pi.

This script checks the construction and computes, for each n and k, I(pi_k),
the exact d(pi_k) (from the certified distance tables, n <= NMAX_EXACT) and the
predicted forced-rotation lower bound n - maxgap, and reports the growth of
d(pi_k) - 2*I(pi_k) - floor(n/2) against n for fixed k.
"""

import argparse
import json
import os
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
sys.path.insert(0, os.path.join(ROOT, "exact"))
from moves import CORE_VERSION  # noqa: E402
from bfs import factorials, rank_perm  # noqa: E402
from line_model import line_bfs  # noqa: E402

VERSION = "h13_blocks-1.1"
OUT = os.path.join(ROOT, "data", "runs", "h13_blocks")


def make_pi(n, k):
    """Identity except k disjoint local swaps, edges at floor(i*n/k), i=0..k-1
    (edge 0 always included). Requires consecutive edges >= 2 apart (disjoint)."""
    edges = sorted(set((i * n) // k for i in range(k)))
    if len(edges) != k:
        return None, None
    for a, b in zip(edges, edges[1:] + [edges[0] + n]):
        if b - a < 2:
            return None, None
    p = list(range(n))
    for e in edges:
        p[e], p[(e + 1) % n] = p[(e + 1) % n], p[e]
    return tuple(p), edges


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def min_inv_line(pi):
    """I(pi) = min over cut q and shift c of inversions of the relabelled line
    (same definition as experiments/line_model.py)."""
    n = len(pi)
    best = None
    for q in range(n):
        for c in range(n):
            w = tuple((pi[(q + 1 + j) % n] - (q + 1 - c)) % n for j in range(n))
            inv = inversions(w)
            if best is None or inv < best:
                best = inv
    return best


def forced_rotation_lb(n, edges):
    """n - maxgap among the k edges placed cyclically on Z_n (0 is one of them)."""
    edges = sorted(edges)
    gaps = [(edges[(i + 1) % len(edges)] - edges[i]) % n or n for i in range(len(edges))]
    return n - max(gaps)


def exact_d(n, pi, fact):
    path = os.path.join(ROOT, "data", "tables", f"dist_n{n}.bin")
    if not os.path.exists(path):
        return None
    with open(path, "rb") as f:
        data = f.read()
    return data[rank_perm(pi, fact)]


def exact_d_line(n, pi, fact):
    """min over (q, c) of the line-restricted distance (as line_model.py)."""
    tables = [line_bfs(n, j1, fact) for j1 in range(n)]
    best = None
    for q in range(n):
        j0 = (-q - 1) % n
        for c in range(n):
            w = tuple((pi[(q + 1 + j) % n] - (q + 1 - c)) % n for j in range(n))
            j1 = (c - q - 1) % n
            d = tables[j1][rank_perm(w, fact) * n + j0]
            if best is None or d < best:
                best = d
    return best


def run(ks, ns, log, args_nmax_dline=10):
    rows = []
    for k in ks:
        for n in ns:
            pi, edges = make_pi(n, k)
            if pi is None:
                log(f"k={k} n={n}: edges not disjoint, skipped")
                continue
            fact = factorials(n)
            I = min_inv_line(pi)
            d = exact_d(n, pi, fact)
            dl = exact_d_line(n, pi, fact) if n <= args_nmax_dline else None
            lb = forced_rotation_lb(n, edges)
            row = {"k": k, "n": n, "edges": edges, "I": I, "d_exact": d, "d_line_exact": dl,
                   "forced_rotation_lb": lb, "floor_n_2": n // 2,
                   "predicted_d_line": I + lb,
                   "gap_over_floor_n2": (None if d is None else d - 2 * I - n // 2),
                   "gap_line_over_floor_n2": (None if dl is None else dl - 2 * I - n // 2)}
            rows.append(row)
            log(f"k={k} n={n}: I(pi)={I} edges={edges} forced_rotation_lb={lb} "
                f"predicted_d_line=I+lb={I + lb} d_line_exact={dl} d_exact={d} "
                f"floor(n/2)={n // 2} d-2I-floor(n/2)={row['gap_over_floor_n2']} "
                f"d_line-2I-floor(n/2)={row['gap_line_over_floor_n2']}")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ks", type=int, nargs="+", default=[2, 3, 4])
    ap.add_argument("--ns", type=int, nargs="+", default=[9, 10, 11, 12])
    ap.add_argument("--nmax_dline", type=int, default=10, help="line BFS is O(n! * n); keep small")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    logf = open(os.path.join(OUT, "report.log"), "a")

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + "\n")
        logf.flush()

    log(f"== {VERSION} {CORE_VERSION} args={vars(args)}")
    rows = run(args.ks, args.ns, log, args.nmax_dline)
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump({"version": VERSION, "core": CORE_VERSION, "rows": rows}, f, indent=1)


if __name__ == "__main__":
    main()
