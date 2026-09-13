"""First step of PLAN §8 candidate (a): trace an actual geodesic (from a
certified distance table, 4 <= n <= 10) from pi down to id_n and record the
sequence of head positions at which X (swap) is used.

Motivation: H10/C27 is REFUTED because the N1 variant (sorting one cycle of
f_c to completion at its carrier, then moving on) loses Theta(n) against B_n
on affine inputs pi(i) = a*i + b (H11, data/runs/h10_verdict/report.md from
PR #5). This script asks a narrower, checkable question: do true geodesics on
these same inputs revisit a head position for more than one X, i.e. do they
interleave work across what the N1 construction treats as separate carriers?
This is diagnostic only -- one arbitrary geodesic per input (greedy first
available decreasing move in a fixed L,R,X tie-break order), not the set of
all geodesics and not a claim about every geodesic.

Usage: python3 experiments/geodesic_trace.py [--n N --a A --b B]...
Requires data/tables/dist_n{n}.bin (rule 9 core tables, n <= 12).
"""

import argparse
import json
import os
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
sys.path.insert(0, os.path.join(ROOT, "exact"))

from moves import apply_move, identity, CORE_VERSION  # noqa: E402
from bfs import factorials, rank_perm  # noqa: E402


def load_table(n):
    fact = factorials(n)
    with open(os.path.join(ROOT, "data", "tables", f"dist_n{n}.bin"), "rb") as f:
        dist = f.read()
    return fact, dist


def geodesic(pi, fact, dist, n):
    """One geodesic id_n -> ... -> pi (built backwards from pi), fixed
    tie-break order X, L, R. Every step is checked to strictly decrease the
    table distance, so the result is a geodesic regardless of tie-break."""
    idn = identity(n)

    def d(p):
        return dist[rank_perm(list(p), fact)]

    cur = pi
    word = []
    while cur != idn:
        dc = d(cur)
        for m in ("X", "L", "R"):
            nxt = apply_move(cur, m)
            if d(nxt) == dc - 1:
                word.append(m)
                cur = nxt
                break
        else:
            raise AssertionError(("no decreasing move", cur))
    return "".join(word)


def trace(word, n):
    """Head position (mod n) at each X in the word, in order; L: head += 1,
    R: head -= 1 (PROBLEM.md conventions; only relative positions matter
    here, so the sign convention does not affect the revisit count)."""
    head = 0
    xs = []
    for ch in word:
        if ch == "L":
            head = (head + 1) % n
        elif ch == "R":
            head = (head - 1) % n
        else:
            xs.append(head)
    return xs


def analyze(pi, n, fact, dist, label):
    Bn = n * (n - 1) // 2
    w = geodesic(pi, fact, dist, n)
    xs = trace(w, n)
    uniq = list(dict.fromkeys(xs))
    revisits = len(xs) - len(uniq)
    return {
        "label": label, "n": n, "pi": list(pi), "d": len(w), "B_n": Bn,
        "word": w, "N_X": w.count("X"), "N_rot": len(w) - w.count("X"),
        "x_head_trace": xs, "distinct_x_heads": len(uniq),
        "head_revisits": revisits,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=10)
    a = ap.parse_args()
    cases = []
    # affine inputs pi(i) = a*i + b (mod n), a != 1: the family that
    # refutes H10 (data/runs/h10_verdict/report.md). One case per n plus the
    # exact C27 part-2 counterexample at n = 10.
    affine = [(9, 2, 0), (10, 7, 2), (10, 3, 4)]
    for n, aa, bb in affine:
        if n > a.nmax:
            continue
        pi = tuple((aa * i + bb) % n for i in range(n))
        fact, dist = load_table(n)
        cases.append(analyze(pi, n, fact, dist, f"affine n={n} (a,b)=({aa},{bb})"))
    part2_n10 = (9, 2, 5, 0, 1, 8, 7, 6, 3, 4)
    if a.nmax >= 10:
        fact, dist = load_table(10)
        cases.append(analyze(part2_n10, 10, fact, dist, "C27 part-2 counterexample n=10"))
    out_dir = os.path.join(ROOT, "data", "runs", "geodesic_trace")
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "report.json"), "w") as f:
        json.dump({"core": CORE_VERSION, "cases": cases}, f, indent=1)
    for c in cases:
        print(f"{c['label']}: d={c['d']} B_n={c['B_n']} N_X={c['N_X']} "
              f"distinct_x_heads={c['distinct_x_heads']} revisits={c['head_revisits']} "
              f"trace={c['x_head_trace']}")


if __name__ == "__main__":
    main()
