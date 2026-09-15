"""toric_cuts-1.0 -- the double-cut grid of a permutation and the candidates
tried for H13-I (session 9).

Model (docs/notes/h13_line_model.md §0, restated in docs/notes/h13i_verdict.md):
for pi in S_n and a cut (q, c) the line is w_j = pi(q+1+j) - (q+1-c) mod n.
Writing a = q + 1 and b = q + 1 - c (a bijection of Z_n^2) this is

    w^{a,b}_j = (pi(a + j) - b) mod n,      F(a, b) = inv(w^{a,b}),
    I(pi) = min_{a,b} F(a, b),              M_n = floor((n-1)^2 / 4).

H13-I: I(pi) <= M_n.

What this module provides:
  cut_grid            -- F by the definition (build every line, count inversions)
  cut_grid_formula    -- F by the closed formula (Lemma 2 of the verdict note)
  mu                  -- the common row/column mean of F (Lemma 3)
  psi                 -- the cyclic correlation of Lemma 4
  circular_swap_dist  -- cs(pi), minimum number of swaps of circularly adjacent
                         positions taking pi to a rotation of the identity
  candidate_report    -- every candidate of the session with its smallest
                         counterexample

Usage:
    python3 experiments/toric_cuts.py --nmax 8 [--cs-nmax 8] [--out DIR]
"""

import argparse
import itertools
import json
import os
import time
from collections import deque
from fractions import Fraction

VERSION = "toric_cuts-1.0"


# ---------------------------------------------------------------- basics

def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def line(perm, a, b):
    """The line of the double cut (a, b)."""
    n = len(perm)
    return tuple((perm[(a + j) % n] - b) % n for j in range(n))


def cut_grid(perm):
    """F(a, b) for all (a, b), straight from the definition."""
    n = len(perm)
    return [[inversions(line(perm, a, b)) for b in range(n)] for a in range(n)]


def cut_grid_formula(perm):
    """F(a, b) by the closed formula of Lemma 2 (verdict note)."""
    n = len(perm)
    ip = [0] * n
    for i, v in enumerate(perm):
        ip[v] = i
    I0 = inversions(perm)
    A = [0] * (n + 1)
    for a in range(n):
        A[a + 1] = A[a] + perm[a]
    B = [0] * (n + 1)
    for b in range(n):
        B[b + 1] = B[b] + ip[b]
    N = [[0] * (n + 1) for _ in range(n + 1)]
    for a in range(n):
        for b in range(n + 1):
            N[a + 1][b] = N[a][b] + (1 if perm[a] < b else 0)
    return [[I0 + (n - 1) * (a + b) + 2 * a * b
             - 2 * A[a] - 2 * B[b] - 2 * n * N[a][b]
             for b in range(n)] for a in range(n)]


def big_T(perm):
    """T(pi) = sum over ordered pairs of ((j-i) mod n)((pi(j)-pi(i)) mod n)."""
    n = len(perm)
    return sum(((j - i) % n) * ((perm[j] - perm[i]) % n)
               for i in range(n) for j in range(n) if i != j)


def mu(perm):
    """Common row/column mean of F, as an exact Fraction (Lemma 3)."""
    n = len(perm)
    return Fraction(n * (n - 1), 2) - Fraction(big_T(perm), n * n)


def psi(perm, a, b):
    """Psi(a, b) = sum_j ((j - a) mod n) * ((pi(j) - b) mod n)  (Lemma 4)."""
    n = len(perm)
    return sum(((j - a) % n) * ((perm[j] - b) % n) for j in range(n))


def I_of(perm):
    return min(min(row) for row in cut_grid_formula(perm))


# ------------------------------------------------- circular swap distance

def circular_swap_table(n):
    """BFS distance to the set of rotations of the identity, using the n swaps
    of circularly adjacent positions.  Returns (perms, index, dist, prev)."""
    perms = list(itertools.permutations(range(n)))
    idx = {p: i for i, p in enumerate(perms)}
    dist = [-1] * len(perms)
    prev = [None] * len(perms)
    dq = deque()
    for c in range(n):
        p = tuple((j + c) % n for j in range(n))
        dist[idx[p]] = 0
        dq.append(p)
    while dq:
        p = dq.popleft()
        d = dist[idx[p]]
        for i in range(n):
            j = (i + 1) % n
            q = list(p)
            q[i], q[j] = q[j], q[i]
            q = tuple(q)
            k = idx[q]
            if dist[k] < 0:
                dist[k] = d + 1
                prev[k] = (p, i)
                dq.append(q)
    return perms, idx, dist, prev


def circular_swap_word(perm, table=None):
    """A shortest list of edges (i means: swap positions i and i+1 mod n)
    sorting perm into a rotation of the identity."""
    n = len(perm)
    perms, idx, dist, prev = table if table else circular_swap_table(n)
    cur = tuple(perm)
    edges = []
    while dist[idx[cur]] > 0:
        p, i = prev[idx[cur]]
        edges.append(i)
        cur = p
    return list(reversed(edges))


# ----------------------------------------------------------- candidates

def representatives(n):
    """One per position-rotation orbit; every toric class is met."""
    for perm in itertools.permutations(range(n)):
        if perm[0] == 0:
            yield perm


def candidate_values(n, perm, G):
    """Each entry: (value, bound) -- the candidate fails when value > bound."""
    M = ((n - 1) ** 2) // 4
    rowmin = [min(G[a][b] for a in range(n)) for b in range(n)]
    out = {
        "H13I": (min(rowmin), M),
        "C37_avg_rowmin": (sum(rowmin), n * M),
        "rowmax": (max(rowmin), M),
        "diag": (min(Fraction(sum(G[a][(a + t) % n] for a in range(n)), n)
                     for t in range(n)), M),
        "anti": (min(Fraction(sum(G[a][(s - a) % n] for a in range(n)), n)
                     for s in range(n)), M),
        "linear": (min(Fraction(sum(G[a][(k * a + t) % n] for a in range(n)), n)
                       for k in range(n) for t in range(n)), M),
        "T1_shift_pos": (min(sum(G[(i + c) % n][perm[i]] for i in range(n))
                             for c in range(n)), n * M),
        "T2_shift_val": (min(sum(G[i][(perm[i] + c) % n] for i in range(n))
                             for c in range(n)), n * M),
        "T3_shift_both": (min(sum(G[(i + c) % n][(perm[i] + d) % n]
                                  for i in range(n))
                              for c in range(n) for d in range(n)), n * M),
    }
    if n % 2 == 0:
        out["T4_antipodal_b"] = (
            max(rowmin[b] + rowmin[(b + n // 2) % n] for b in range(n)), 2 * M)
    return out


def candidate_report(nmax):
    """Smallest n and lexicographically first representative refuting each
    candidate; candidates with no counterexample up to nmax are listed too."""
    first = {}
    checked = {}
    for n in range(4, nmax + 1):
        t0 = time.time()
        cnt = 0
        for perm in representatives(n):
            cnt += 1
            G = cut_grid_formula(perm)
            for key, (val, bound) in candidate_values(n, perm, G).items():
                checked.setdefault(key, set()).add(n)
                if val > bound and key not in first:
                    first[key] = {"n": n, "perm": list(perm),
                                  "value": str(val), "bound": bound}
        print("  n=%d: %d representatives, %.1f s" % (n, cnt, time.time() - t0))
    survivors = sorted(k for k in checked if k not in first)
    return {"refuted": first, "not_refuted_up_to_nmax": survivors,
            "nmax": nmax}


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--cs-nmax", type=int, default=8)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    out = args.out or os.path.join(root, "data", "runs", "h13i_verdict")
    os.makedirs(out, exist_ok=True)

    print("%s: candidate scan up to n=%d" % (VERSION, args.nmax))
    rep = candidate_report(args.nmax)

    print("circular swap distance vs I, up to n=%d" % args.cs_nmax)
    cs = []
    for n in range(4, args.cs_nmax + 1):
        t0 = time.time()
        perms, idx, dist, prev = circular_swap_table(n)
        M = ((n - 1) ** 2) // 4
        maxcs = 0
        diff = 0
        first = None      # lexicographically first pi with cs < I
        maxgap = 0
        for p in perms:
            I = I_of(p)
            c = dist[idx[p]]
            assert c <= I, (n, p, c, I)
            maxcs = max(maxcs, c)
            if c != I:
                diff += 1
                maxgap = max(maxgap, I - c)
                if first is None:
                    first = {"perm": list(p), "cs": c, "I": I, "gap": I - c,
                             "swap_edges": circular_swap_word(
                                 p, (perms, idx, dist, prev))}
        cs.append({"n": n, "M": M, "max_cs": maxcs, "num_cs_lt_I": diff,
                   "max_gap": maxgap, "lex_first_cs_lt_I": first,
                   "seconds": round(time.time() - t0, 2)})
        print("  n=%d max cs=%d (M=%d), cs<I on %d of %d, max gap=%d, first=%s"
              % (n, maxcs, M, diff, len(perms), maxgap, first))

    res = {"version": VERSION, "candidates": rep, "circular_swap": cs}
    with open(os.path.join(out, "toric_cuts.json"), "w") as f:
        json.dump(res, f, indent=1)
    print("written:", os.path.join(out, "toric_cuts.json"))


if __name__ == "__main__":
    main()
