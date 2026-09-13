"""Rotation budget in the line model (H13, session 8).

Line model (experiments/line_model.py): circle positions q+1..q+n form a line,
the head never crosses or swaps across the cut edge {q, q+1}; with shift c the
relabelled line w must become the identity, the head starts at line index
j0 = (-q-1) mod n and must end at j1 = (c-q-1) mod n.

A *reduced* line word swaps only inverted adjacent pairs, so N_X = inv(w).
Quantities (all exhaustive over pi, 4 <= n <= NMAX):
  R(w, j0, j1)  = min N_rot over reduced line words for (w, j0, j1)
                  (0-1 BFS backwards from (id, j1): rotations cost 1, swaps 0);
  d_red(pi)     = min over (q, c) of inv(w) + R(w, j0, j1)   (>= d_line >= d);
  Rbest(pi)     = min over (q, c) with inv(w) = I(pi) of R(w, j0, j1)
                  (rotation budget when the swap count is the minimum I(pi));
  Rfree(w)      = min over j0, j1 of R(w, j0, j1)  (line only, no cut choice).
Reported: max_pi d_red - d, max_pi d_red - I - floor(n^2/4), max_pi Rbest,
max_w Rfree, histograms, argmax examples.
Usage: python3 experiments/line_rotation_budget.py --nmax 8
Output: data/runs/line_rotation_budget/.  Version line_rotation_budget-1.0.
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

VERSION = "line_rotation_budget-1.0"
OUT = os.path.join(ROOT, "data", "runs", "line_rotation_budget")


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def reduced_line_01bfs(n, j1, fact):
    """R(w, j0, j1) for all (w, j0): 0-1 BFS backwards from (id, j1).
    Backward move = inverse of a forward move: a forward swap at edge p
    (allowed iff w[p] > w[p+1]) is undone by swapping at p when w[p] < w[p+1]
    (it creates an inversion in the backward direction).  Rotations cost 1."""
    total = fact[n]
    INF = 255
    dist = bytearray([INF]) * (total * n)
    start = (tuple(range(n)), j1)
    dist[rank_perm(start[0], fact) * n + j1] = 0
    dq = deque([start])
    while dq:
        w, p = dq.popleft()
        idx = rank_perm(w, fact) * n + p
        d = dist[idx]
        # swap (cost 0): backwards allowed iff w[p] < w[p+1]
        if p < n - 1 and w[p] < w[p + 1]:
            w2 = list(w)
            w2[p], w2[p + 1] = w2[p + 1], w2[p]
            w2 = tuple(w2)
            i2 = rank_perm(w2, fact) * n + p
            if dist[i2] > d:
                dist[i2] = d
                dq.appendleft((w2, p))
        for p2 in (p - 1, p + 1):
            if 0 <= p2 < n:
                i2 = rank_perm(w, fact) * n + p2
                if dist[i2] > d + 1:
                    dist[i2] = d + 1
                    dq.append((w, p2))
    return dist


def run(n, log, line_tables=None):
    t0 = time.time()
    fact = factorials(n)
    tabs = [reduced_line_01bfs(n, j1, fact) for j1 in range(n)]
    log(f"n={n}: {n} reduced 0-1 BFS runs, {time.time() - t0:.0f} s")
    path = os.path.join(ROOT, "data", "tables", f"dist_n{n}.bin")
    circ = open(path, "rb").read()
    B = n * (n - 1) // 2
    q2 = n * n // 4
    # Rfree over lines
    max_rfree, arg_rfree = -1, None
    hist_rfree = {}
    for w in itertools.permutations(range(n)):
        r = rank_perm(w, fact)
        best = min(tabs[j1][r * n + j0] for j1 in range(n) for j0 in range(n))
        hist_rfree[best] = hist_rfree.get(best, 0) + 1
        if best > max_rfree:
            max_rfree, arg_rfree = best, w
    log(f"n={n}: max Rfree over lines = {max_rfree} at {arg_rfree}; floor(n^2/4) = {q2}")
    max_gap, arg_gap = -1, None
    hist_gap = {}
    max_ex, arg_ex = -99, None
    max_rbest, arg_rbest = -1, None
    hist_rbest = {}
    hist_ex = {}
    prof = {}
    for pi in itertools.permutations(range(n)):
        dc = circ[rank_perm(pi, fact)]
        best_len, I, rbest = None, None, None
        cand = []
        for q in range(n):
            j0 = (-q - 1) % n
            for c in range(n):
                w = tuple((pi[(q + 1 + j) % n] - (q + 1 - c)) % n for j in range(n))
                j1 = (c - q - 1) % n
                inv = inversions(w)
                rot = tabs[j1][rank_perm(w, fact) * n + j0]
                cand.append((inv, rot))
        I = min(inv for inv, rot in cand)
        rbest = min(rot for inv, rot in cand if inv == I)
        best_len = min(inv + rot for inv, rot in cand)
        gap = best_len - dc
        hist_gap[gap] = hist_gap.get(gap, 0) + 1
        if gap > max_gap:
            max_gap, arg_gap = gap, pi
        ex = best_len - I - q2
        hist_ex[ex] = hist_ex.get(ex, 0) + 1
        if ex > max_ex:
            max_ex, arg_ex = ex, pi
        hist_rbest[rbest] = hist_rbest.get(rbest, 0) + 1
        if rbest > max_rbest:
            max_rbest, arg_rbest = rbest, pi
        p = prof.setdefault(I, {"count": 0, "max_d_red": 0, "max_rbest": 0, "max_d": 0})
        p["count"] += 1
        p["max_d_red"] = max(p["max_d_red"], best_len)
        p["max_rbest"] = max(p["max_rbest"], rbest)
        p["max_d"] = max(p["max_d"], dc)
    row = {"n": n, "B_n": B, "floor_n2_4": q2, "seconds": round(time.time() - t0),
           "max_Rfree_line": max_rfree, "argmax_Rfree_line": list(arg_rfree),
           "hist_Rfree_line": dict(sorted(hist_rfree.items())),
           "max_d_red_minus_d": max_gap, "argmax_d_red_minus_d": list(arg_gap),
           "hist_d_red_minus_d": dict(sorted(hist_gap.items())),
           "max_d_red_minus_I_minus_floor_n2_4": max_ex, "argmax": list(arg_ex),
           "hist_d_red_minus_I_minus_floor_n2_4": dict(sorted(hist_ex.items())),
           "max_Rbest": max_rbest, "argmax_Rbest": list(arg_rbest),
           "hist_Rbest": dict(sorted(hist_rbest.items())),
           "profile_by_I": {k: prof[k] for k in sorted(prof)}}
    log(f"n={n}: max d_red - d = {max_gap} at {arg_gap}, hist {dict(sorted(hist_gap.items()))}")
    log(f"n={n}: max d_red - I - floor(n^2/4) = {max_ex} at {arg_ex}, hist {dict(sorted(hist_ex.items()))}")
    log(f"n={n}: max Rbest = {max_rbest} at {arg_rbest}, hist {dict(sorted(hist_rbest.items()))}")
    for I in sorted(prof):
        p = prof[I]
        log(f"   I={I:2d} count={p['count']:6d} max_d={p['max_d']:2d} max_d_red={p['max_d_red']:2d} "
            f"max_Rbest={p['max_rbest']:2d} I+floor(n^2/4)={I + q2}")
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
        with open(os.path.join(OUT, f"report_n{args.nmin}_{n}.json"), "w") as f:
            json.dump({"version": VERSION, "core": CORE_VERSION, "rows": rows}, f, indent=1)


if __name__ == "__main__":
    main()
