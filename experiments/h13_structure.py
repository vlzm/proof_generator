"""H13 step (2), parts (1) and (3): structure of d_line - d and of near-maximal
I(pi), exhaustive at small n (reuses experiments/line_model.py's line BFS).

For every pi at 4 <= n <= NMAX:
  - I(pi) = min_{q,c} inversions of the relabelled line (as in line_model.py);
  - d_line(pi), d(pi) from the certified tables;
  - blocks(pi) = number of maximal runs of non-fixed positions in the best
    (q, c) relabelled line (a coarse measure of how "scattered" the mismatches
    of pi are; the h13_blocks.py family has blocks(pi) = k).

Reports: (1) the gap d_line - d grouped by blocks; (3) the permutations with
I(pi) > floor((n-1)^2/4) - n/4, and their blocks/case (reflection or not).
"""

import argparse
import itertools
import json
import os
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
sys.path.insert(0, os.path.join(ROOT, "exact"))
from moves import CORE_VERSION  # noqa: E402
from bfs import factorials, rank_perm  # noqa: E402
from line_model import line_bfs, inversions  # noqa: E402

VERSION = "h13_structure-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h13_structure")


def blocks_of(w):
    """Number of maximal runs of positions j with w[j] != j (fixed-point gaps
    split the line into blocks); w is the relabelled line array."""
    n = len(w)
    cnt = 0
    prev_bad = False
    for j in range(n):
        bad = w[j] != j
        if bad and not prev_bad:
            cnt += 1
        prev_bad = bad
    return cnt


def run(n, log):
    fact = factorials(n)
    tables = [line_bfs(n, j1, fact) for j1 in range(n)]
    path = os.path.join(ROOT, "data", "tables", f"dist_n{n}.bin")
    circ = open(path, "rb").read()
    B = n * (n - 1) // 2
    thresh = (n - 1) ** 2 // 4 - n / 4
    gap_by_blocks = {}
    near_max = []
    for pi in itertools.permutations(range(n)):
        best_dl, best_inv, best_blocks = None, None, None
        for q in range(n):
            j0 = (-q - 1) % n
            for c in range(n):
                w = tuple((pi[(q + 1 + j) % n] - (q + 1 - c)) % n for j in range(n))
                j1 = (c - q - 1) % n
                d = tables[j1][rank_perm(w, fact) * n + j0]
                if best_dl is None or d < best_dl:
                    best_dl = d
                inv = inversions(w)
                if best_inv is None or inv < best_inv:
                    best_inv = inv
                    best_blocks = blocks_of(w)
        dc = circ[rank_perm(pi, fact)]
        gap = best_dl - dc
        gap_by_blocks.setdefault(best_blocks, {}).setdefault(gap, 0)
        gap_by_blocks[best_blocks][gap] += 1
        if best_inv > thresh:
            near_max.append({"pi": list(pi), "I": best_inv, "blocks": best_blocks,
                              "d_line": best_dl, "d": int(dc)})
    log(f"n={n}: floor((n-1)^2/4)={(n - 1) ** 2 // 4}, threshold(>)={thresh}, "
        f"near_max count={len(near_max)}")
    log(f"n={n}: gap(d_line-d) by blocks: "
        f"{ {k: dict(sorted(v.items())) for k, v in sorted(gap_by_blocks.items())} }")
    return {"n": n, "B_n": B, "threshold": thresh, "near_max": near_max,
            "gap_by_blocks": {k: v for k, v in gap_by_blocks.items()}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=7)
    ap.add_argument("--nmax", type=int, default=7)
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
        with open(os.path.join(OUT, f"report_n{args.nmin}_{n}.json"), "w") as f:
            json.dump({"version": VERSION, "core": CORE_VERSION, "rows": rows}, f, indent=1)


if __name__ == "__main__":
    main()
