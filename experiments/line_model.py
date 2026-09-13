"""Line model (H13 candidate): does one unused edge ("cut") of the circle suffice?

For a permutation pi, a cut edge {q, q+1} and a shift c, sort pi on the LINE of
circle positions q+1, ..., q+n (the head never crosses the cut and never swaps
across it).  d_line(pi) = min over (q, c) of the length of a shortest such word.
Trivially d_line(pi) >= d(pi).  On reflections the two-arc word of C32 uses one
cut only, so d_line(pi_h) = d(pi_h) = B_n - delta_n(h, 1).

Question: is max_pi d_line(pi) = B_n (i.e. is the wrap-around never needed for
the diameter), and how often d_line(pi) = d(pi)?  Also the purely combinatorial
part: I(pi) = min over (q, c) of the number of inversions of the relabelled line
(a lower bound on N_X in the line model); is max_pi I(pi) = floor((n-1)^2/4)?

Method: BFS in the line graph (states: arrangement of n elements on the line,
head index 0..n-1; moves: head +-1 inside the line, swap at head < n-1) from
(identity, head j1) for every j1; then for every pi, q, c look up the distance
from the relabelled start state.  Exhaustive for 4 <= n <= NMAX.
Usage: python3 experiments/line_model.py --nmax 8
Output: data/runs/line_model/report.json, report.md.  Version line_model-1.0.

Version 1.1 (session 8) adds the intermediate "walk" model: the head may cross
the cut (moves +-1 mod n) but never swaps across it (no swap at line index
n-1).  d <= d_walk <= d_line, so d_line - d splits into the part saved by
walking around the circle (d_line - d_walk) and the part saved by swapping
across the cut (d_walk - d).  With --walk the per-pi joint histogram of
(d_line - d_walk, d_walk - d) is reported, also grouped by I(pi) (read from
data/runs/line_profile/I_n{n}.bin if present) and by d.
"""

import argparse
import itertools
import json
import os
import sys
import time
from collections import deque

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
sys.path.insert(0, os.path.join(ROOT, "exact"))
from moves import CORE_VERSION  # noqa: E402
from bfs import factorials, rank_perm  # noqa: E402

VERSION = "line_model-1.1"
OUT = os.path.join(ROOT, "data", "runs", "line_model")


def line_bfs(n, j1, fact):
    """Distances in the line graph from (identity, head j1); index = rank*n + head."""
    N = fact[n] * n if len(fact) > n else None
    total = 1
    for k in range(2, n + 1):
        total *= k
    dist = bytearray([255]) * (total * n)
    start = (tuple(range(n)), j1)
    dist[rank_perm(start[0], fact) * n + j1] = 0
    dq = deque([start])
    while dq:
        w, p = dq.popleft()
        d = dist[rank_perm(w, fact) * n + p]
        nxt = []
        if p > 0:
            nxt.append((w, p - 1))
        if p < n - 1:
            nxt.append((w, p + 1))
            w2 = list(w)
            w2[p], w2[p + 1] = w2[p + 1], w2[p]
            nxt.append((tuple(w2), p))
        for w2, p2 in nxt:
            idx = rank_perm(w2, fact) * n + p2
            if dist[idx] == 255:
                dist[idx] = d + 1
                dq.append((w2, p2))
    return dist


def walk_bfs(n, j1, fact):
    """Distances in the walk graph (head moves +-1 mod n, swap only at head < n-1)
    from (identity, head j1); index = rank*n + head."""
    total = 1
    for k in range(2, n + 1):
        total *= k
    dist = bytearray([255]) * (total * n)
    start = (tuple(range(n)), j1)
    dist[rank_perm(start[0], fact) * n + j1] = 0
    dq = deque([start])
    while dq:
        w, p = dq.popleft()
        d = dist[rank_perm(w, fact) * n + p]
        nxt = [(w, (p - 1) % n), (w, (p + 1) % n)]
        if p < n - 1:
            w2 = list(w)
            w2[p], w2[p + 1] = w2[p + 1], w2[p]
            nxt.append((tuple(w2), p))
        for w2, p2 in nxt:
            idx = rank_perm(w2, fact) * n + p2
            if dist[idx] == 255:
                dist[idx] = d + 1
                dq.append((w2, p2))
    return dist


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def run(n, log, walk=False):
    t0 = time.time()
    fact = factorials(n)
    tables = [line_bfs(n, j1, fact) for j1 in range(n)]
    log(f"n={n}: {n} line BFS runs, {time.time() - t0:.0f} s")
    wtables = None
    if walk:
        wtables = [walk_bfs(n, j1, fact) for j1 in range(n)]
        log(f"n={n}: {n} walk BFS runs, {time.time() - t0:.0f} s")
    ipath = os.path.join(ROOT, "data", "runs", "line_profile", f"I_n{n}.bin")
    Ivals = open(ipath, "rb").read() if os.path.exists(ipath) else None
    joint = {}          # (d_line - d_walk, d_walk - d) -> count
    by_I = {}           # I -> {(gl, gw): count}
    by_d = {}           # d -> {(gl, gw): count}
    walk_max, walk_arg = 0, None
    B = n * (n - 1) // 2
    path = os.path.join(ROOT, "data", "tables", f"dist_n{n}.bin")
    circ = open(path, "rb").read() if os.path.exists(path) else None
    max_dl, argmax_dl, eq_count, hist_gap = 0, None, 0, {}
    max_inv, argmax_inv = 0, None
    cnt = 0
    for pi in itertools.permutations(range(n)):
        cnt += 1
        best, best_inv, best_walk = None, None, None
        for q in range(n):
            j0 = (-q - 1) % n
            for c in range(n):
                # line index j <-> circle position q+1+j; target line[j] = q+1+j-c
                # relabel element e -> e - (q+1-c) so that the target is identity
                w = tuple((pi[(q + 1 + j) % n] - (q + 1 - c)) % n for j in range(n))
                j1 = (c - q - 1) % n
                d = tables[j1][rank_perm(w, fact) * n + j0]
                if best is None or d < best:
                    best = d
                if wtables is not None:
                    dw = wtables[j1][rank_perm(w, fact) * n + j0]
                    if best_walk is None or dw < best_walk:
                        best_walk = dw
                inv = inversions(w)
                if best_inv is None or inv < best_inv:
                    best_inv = inv
        if best > max_dl:
            max_dl, argmax_dl = best, pi
        if best_inv > max_inv:
            max_inv, argmax_inv = best_inv, pi
        if circ is not None:
            rk = rank_perm(pi, fact)
            dc = circ[rk]
            gap = best - dc
            hist_gap[gap] = hist_gap.get(gap, 0) + 1
            if gap == 0:
                eq_count += 1
            if wtables is not None:
                key = (best - best_walk, best_walk - dc)
                joint[key] = joint.get(key, 0) + 1
                if best_walk - dc > walk_max:
                    walk_max, walk_arg = best_walk - dc, pi
                if Ivals is not None:
                    h = by_I.setdefault(Ivals[rk], {})
                    h[key] = h.get(key, 0) + 1
                h = by_d.setdefault(dc, {})
                h[key] = h.get(key, 0) + 1
    row = {"n": n, "B_n": B, "count": cnt, "max_d_line": max_dl, "argmax_d_line": list(argmax_dl),
           "max_min_inversions": max_inv, "floor_(n-1)^2/4": (n - 1) ** 2 // 4,
           "argmax_inv": list(argmax_inv), "d_line_eq_d": eq_count,
           "hist_d_line_minus_d": dict(sorted(hist_gap.items())), "seconds": round(time.time() - t0)}
    if wtables is not None:
        fmt = lambda h: {f"{a},{b}": v for (a, b), v in sorted(h.items())}  # noqa: E731
        row["walk"] = {"joint_hist_(d_line-d_walk, d_walk-d)": fmt(joint),
                       "max_d_walk_minus_d": walk_max, "argmax_d_walk_minus_d": list(walk_arg) if walk_arg else None,
                       "by_I": {int(k): fmt(v) for k, v in sorted(by_I.items())},
                       "by_d": {int(k): fmt(v) for k, v in sorted(by_d.items())}}
        log(f"n={n}: joint hist (d_line - d_walk, d_walk - d): {fmt(joint)}; max d_walk - d = {walk_max} at {walk_arg}")
        for k in sorted(by_I):
            log(f"   I={k:2d}: {fmt(by_I[k])}")
    log(f"n={n}: max d_line = {max_dl} (B_n = {B}) at {argmax_dl}; max min-inv = {max_inv} "
        f"(floor((n-1)^2/4) = {(n - 1) ** 2 // 4}) at {argmax_inv}; d_line = d on {eq_count}/{cnt}; "
        f"hist d_line - d: {dict(sorted(hist_gap.items()))}; {time.time() - t0:.0f} s")
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=4)
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--walk", action="store_true", help="also compute the walk model (head may cross the cut)")
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
        rows.append(run(n, log, walk=args.walk))
        suffix = "_walk" if args.walk else ""
        with open(os.path.join(OUT, f"report_n{args.nmin}_{n}{suffix}.json"), "w") as f:
            json.dump({"version": VERSION, "core": CORE_VERSION, "rows": rows}, f, indent=1)


if __name__ == "__main__":
    main()
