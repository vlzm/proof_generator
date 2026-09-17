"""circular_sorting.py — circular sorting time cs(pi) versus I(pi) (session 9).

cs(pi) — the minimum number of cyclically adjacent transpositions (swap the
contents of positions e and e+1 mod n, any e, no head) taking pi to a rotation
of the identity.  This is the quantity of the paper

  R. M. Adin, N. Alon, Y. Roichman, "Circular sorting", arXiv:2502.14398 (2025),

whose reported main theorem is max_pi cs(pi) = floor((n-1)^2/4) (status CLAIMED
in CLAIMS.md: the full text was not accessible in session 9).

Facts checked here by exhaustive BFS over S_n (start set: the n rotations of the
identity; generators: the n cyclically adjacent transpositions):

  * max_pi cs(pi) = floor((n-1)^2/4)  (the same value as max_pi I(pi), H13-I);
  * cs(pi) <= I(pi) for every pi (bubble sort along the best double cut never
    uses the cut edge), with equality for n <= 7 but not for n = 8, 9: sorting
    that uses all n edges can beat every cut.

cs(pi) is a lower bound for N_X of any word of the LRX model that sorts pi.

Usage: python3 experiments/circular_sorting.py [--nmax 9]
Output: data/runs/circular_sorting/circular_sorting.json.
Version circular_sorting-1.0.
"""

import argparse
import itertools
import json
import os
import time
from collections import deque

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
VERSION = "circular_sorting-1.0"


def inv_table(pi):
    """inv(a, s) for all n^2 double cuts (C37 Lemma 7)."""
    n = len(pi)
    pinv = [0] * n
    for i, v in enumerate(pi):
        pinv[v] = i
    cur_a = sum(1 for i in range(n) for j in range(i + 1, n) if pi[i] > pi[j])
    tab = []
    for a in range(n):
        if a:
            cur_a += n - 1 - 2 * pi[a - 1]
        row, cur = [cur_a], cur_a
        for s in range(1, n):
            cur += n - 1 - 2 * ((pinv[s - 1] - a) % n)
            row.append(cur)
        tab.append(row)
    return tab


def scan(n):
    t0 = time.time()
    perms = list(itertools.permutations(range(n)))
    idx = {p: i for i, p in enumerate(perms)}
    dist = [-1] * len(perms)
    q = deque()
    for r in range(n):
        p = tuple((i + r) % n for i in range(n))
        dist[idx[p]] = 0
        q.append(p)
    while q:
        p = q.popleft()
        d = dist[idx[p]]
        for e in range(n):
            lp = list(p)
            lp[e], lp[(e + 1) % n] = lp[(e + 1) % n], lp[e]
            j = idx[tuple(lp)]
            if dist[j] < 0:
                dist[j] = d + 1
                q.append(tuple(lp))
    max_cs, arg_cs, differ, worst_gap, arg_gap, cs_gt_I = 0, None, 0, 0, None, 0
    for p in perms:
        I = min(min(r) for r in inv_table(p))
        c = dist[idx[p]]
        if c > max_cs:
            max_cs, arg_cs = c, list(p)
        if c > I:
            cs_gt_I += 1
        if c != I:
            differ += 1
            if I - c > worst_gap:
                worst_gap, arg_gap = I - c, list(p)
    return {"n": n, "perms": len(perms), "max_cs": max_cs, "target": (n - 1) ** 2 // 4,
            "argmax_cs": arg_cs, "cs_ne_I": differ, "cs_gt_I": cs_gt_I,
            "worst_I_minus_cs": worst_gap, "argmax_gap": arg_gap,
            "seconds": round(time.time() - t0, 1)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=4)
    ap.add_argument("--nmax", type=int, default=9)
    ap.add_argument("--out", default=os.path.join(ROOT, "data", "runs", "circular_sorting"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rows = []
    for n in range(args.nmin, args.nmax + 1):
        r = scan(n)
        rows.append(r)
        print(json.dumps(r), flush=True)
    path = os.path.join(args.out, "circular_sorting.json")
    old = []
    if os.path.exists(path):
        with open(path) as f:
            old = json.load(f).get("rows", [])
    keep = [r for r in old if not any(r["n"] == x["n"] for x in rows)]
    with open(path, "w") as f:
        json.dump({"version": VERSION, "rows": sorted(keep + rows, key=lambda r: r["n"])}, f, indent=1)
    print("written", path)


if __name__ == "__main__":
    main()
