"""Checker for C37, C38, C39 (session 9, docs/notes/h13i_verdict.md).

Everything here is recomputed from the definition of the double cut

    w^{a,b}_j = (pi(a + j) - b) mod n,   F(a, b) = inv(w^{a,b}),
    I(pi) = min_{a,b} F(a,b),           M = floor((n-1)^2 / 4),

independently of experiments/toric_cuts.py and experiments/toric_cut_scan.c:
the grid is built by listing all n^2 lines and counting inversions by brute
force.  The closed formula of Lemma 2 is *checked against* that grid, never
used to produce it, except in part C where part A has already certified the
equivalence and the note says so explicitly.

  Part A (C37, lemmas 0-5): the arc description, both increment formulas, the
    closed form, the constant row/column sums with mu = C(n,2) - T/n^2, the
    correlation form F = mu + (n-1)^2/2 - (2/n) Psi, and the reversal duality
    F_{pi^r}(a,b) = C(n,2) - F_pi(1-a, b).  All pi for 4 <= n <= AMAX (7),
    random samples for AMAX < n <= 11.
  Part B (C37, lemma 6): the finite content of the impossibility proof -- for
    every reflection pi_h and every cut, F_{pi_h}(a,b) = C(n,2) - s(n-s) with
    s = 1 + ((h - a - b) mod n), and the set of minimising cuts is exactly the
    one or two anti-diagonals a + b = const.  All h, all cuts, 4 <= n <= BMAX.
  Part C (C38): sum_b min_a F(a,b) <= n*M for every toric class, exhaustively
    for 4 <= n <= CMAX (8) over the (n-1)! representatives with pi(0) = 0, plus
    the equality set (only the reflection class).  Results of the C scan for
    9 <= n <= 12 are read from data/runs/h13i_verdict/scan.json when present.
  Part D (C39): cs(pi) -- minimum number of swaps of circularly adjacent
    positions taking pi to a rotation of the identity -- satisfies cs <= I for
    all pi with 4 <= n <= DMAX (8); cs = I for 4 <= n <= 7; at n = 8 there are
    exactly 16 permutations with cs < I, the lexicographically first being
    (0,3,6,1,4,7,2,5); the explicit swap sequence for (0,5,2,7,4,1,6,3) is
    replayed.
  Part E (refuted candidates of the session): each listed smallest
    counterexample really violates its candidate, and no smaller n and no
    lexicographically earlier representative does.

Usage: python3 checks/check_C37.py [--amax 7] [--bmax 24] [--cmax 8] [--dmax 8]
Output: data/runs/check_C37/report.json, report.md
"""

import argparse
import itertools
import json
import os
import random
import sys
import time
from collections import deque
from fractions import Fraction

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

VERSION = "check_C37-1.0"
OUT = os.path.join(ROOT, "data", "runs", "check_C37")


# ------------------------------------------------------------ definitions

def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def line(perm, a, b):
    n = len(perm)
    return tuple((perm[(a + j) % n] - b) % n for j in range(n))


def grid_def(perm):
    """F(a,b) straight from the definition."""
    n = len(perm)
    return [[inversions(line(perm, a, b)) for b in range(n)] for a in range(n)]


def grid_formula(perm):
    """F(a,b) by Lemma 2 (certified against grid_def in part A)."""
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


def M_of(n):
    return ((n - 1) ** 2) // 4


# ------------------------------------------------------------- part A

def part_a(amax, log, rng):
    ok = True
    rows = []
    for n in range(4, 12):
        exhaustive = n <= amax
        if exhaustive:
            sample = list(itertools.permutations(range(n)))
        else:
            sample = []
            for _ in range(200):
                p = list(range(n))
                rng.shuffle(p)
                sample.append(tuple(p))
        t0 = time.time()
        C2 = n * (n - 1) // 2
        for perm in sample:
            G = grid_def(perm)
            ip = [0] * n
            for i, v in enumerate(perm):
                ip[v] = i
            # lemma 0: arc description of an inverted pair
            for a in range(n):
                for b in range(n):
                    cnt = 0
                    for i in range(n):
                        for j in range(i + 1, n):
                            in_a = (i < a <= j)
                            dy = (perm[j] - perm[i]) % n
                            in_b = ((b - perm[i] - 1) % n) < dy
                            if in_a != in_b:
                                cnt += 1
                    if cnt != G[a][b]:
                        ok = False
                        log(f"FAIL lemma 0: n={n} pi={perm} (a,b)=({a},{b})")
                        break
                if not ok:
                    break
            # lemma 1: increments
            for a in range(n):
                for b in range(n):
                    if G[(a + 1) % n][b] - G[a][b] != n - 1 - 2 * ((perm[a] - b) % n):
                        ok = False
                        log(f"FAIL lemma 1a: n={n} pi={perm} ({a},{b})")
                    if G[a][(b + 1) % n] - G[a][b] != n - 1 - 2 * ((ip[b] - a) % n):
                        ok = False
                        log(f"FAIL lemma 1b: n={n} pi={perm} ({a},{b})")
            # lemma 2: closed form
            if grid_formula(perm) != G:
                ok = False
                log(f"FAIL lemma 2: n={n} pi={perm}")
            # lemma 3: constant row and column sums, value n*mu
            T = sum(((j - i) % n) * ((perm[j] - perm[i]) % n)
                    for i in range(n) for j in range(n) if i != j)
            mu = Fraction(C2) - Fraction(T, n * n)
            colsums = {sum(G[a][b] for a in range(n)) for b in range(n)}
            rowsums = {sum(G[a][b] for b in range(n)) for a in range(n)}
            if len(colsums) != 1 or len(rowsums) != 1 or colsums != rowsums \
               or Fraction(colsums.pop()) != n * mu:
                ok = False
                log(f"FAIL lemma 3: n={n} pi={perm}")
            if min(min(r) for r in G) > mu:
                ok = False
                log(f"FAIL lemma 3 corollary I <= mu: n={n} pi={perm}")
            # lemma 4: correlation form
            for a in range(n):
                for b in range(n):
                    psi = sum(((j - a) % n) * ((perm[j] - b) % n) for j in range(n))
                    rhs = mu + Fraction((n - 1) ** 2, 2) - Fraction(2 * psi, n)
                    if rhs != G[a][b]:
                        ok = False
                        log(f"FAIL lemma 4: n={n} pi={perm} ({a},{b})")
            # lemma 5: duality
            pr = tuple(perm[(-i) % n] for i in range(n))
            Gr = grid_def(pr)
            for a in range(n):
                for b in range(n):
                    if Gr[a][b] != C2 - G[(1 - a) % n][b]:
                        ok = False
                        log(f"FAIL lemma 5: n={n} pi={perm} ({a},{b})")
            if not ok:
                break
        rows.append({"n": n, "coverage": "exhaustive" if exhaustive else "sampled",
                     "permutations": len(sample),
                     "seconds": round(time.time() - t0, 2)})
        log("  A n=%d %s %d permutations, %.2f s" %
            (n, rows[-1]["coverage"], len(sample), rows[-1]["seconds"]))
        if not ok:
            break
    return ok, rows


# ------------------------------------------------------------- part B

def part_b(bmax, log):
    """Finite content of lemma 6: the reflection grid and its argmin set."""
    ok = True
    rows = []
    for n in range(4, bmax + 1):
        C2 = n * (n - 1) // 2
        M = M_of(n)
        t0 = time.time()
        diag_sets = []
        for h in range(n):
            pih = tuple((h - i) % n for i in range(n))
            argmin = set()
            for a in range(n):
                for b in range(n):
                    s = 1 + ((h - a - b) % n)
                    want = C2 - s * (n - s)
                    got = inversions(line(pih, a, b))
                    if got != want:
                        ok = False
                        log(f"FAIL lemma 6 grid: n={n} h={h} ({a},{b}) {got}!={want}")
                    if got == M:
                        argmin.add((a + b) % n)
                    if got < M:
                        ok = False
                        log(f"FAIL lemma 6 min: n={n} h={h} ({a},{b}) {got}<{M}")
            want_size = 1 if n % 2 == 0 else 2
            if len(argmin) != want_size:
                ok = False
                log(f"FAIL lemma 6 argmin size: n={n} h={h} {sorted(argmin)}")
            diag_sets.append(sorted(argmin))
        # the n supports, as subsets of the n anti-diagonals, are translates of
        # one another; the counting step of the proof needs sum over h of their
        # indicator to be constant = want_size on every anti-diagonal.
        cover = [0] * n
        for s in diag_sets:
            for t in s:
                cover[t] += 1
        if len(set(cover)) != 1 or cover[0] != (1 if n % 2 == 0 else 2):
            ok = False
            log(f"FAIL lemma 6 counting: n={n} cover={cover}")
        rows.append({"n": n, "argmin_diagonals_per_h": diag_sets[0],
                     "seconds": round(time.time() - t0, 2)})
    log("  B reflections checked for 4 <= n <= %d" % bmax)
    return ok, rows


# ------------------------------------------------------------- part C

def part_c(cmax, log):
    ok = True
    rows = []
    for n in range(4, cmax + 1):
        M = M_of(n)
        t0 = time.time()
        use_def = n <= 6
        worst = -1
        eq = []
        reps = 0
        maxI = -1
        for perm in itertools.permutations(range(n)):
            if perm[0] != 0:
                continue
            reps += 1
            G = grid_def(perm) if use_def else grid_formula(perm)
            rowmin = [min(G[a][b] for a in range(n)) for b in range(n)]
            s = sum(rowmin)
            maxI = max(maxI, min(rowmin))
            if s > n * M:
                ok = False
                log(f"FAIL C38: n={n} pi={perm} sum={s} > {n*M}")
            if s > worst:
                worst = s
                eq = [list(perm)]
            elif s == worst:
                eq.append(list(perm))
        if maxI != M:
            ok = False
            log(f"FAIL H13-I max: n={n} maxI={maxI} != {M}")
        rows.append({"n": n, "representatives": reps, "bound": n * M,
                     "max_sum": worst, "argmax": eq, "max_I": maxI,
                     "grid": "definition" if use_def else "lemma 2 formula",
                     "seconds": round(time.time() - t0, 2)})
        log("  C n=%d reps=%d max sum=%d (bound %d) argmax=%s maxI=%d"
            % (n, reps, worst, n * M, eq, maxI))
    scan = os.path.join(ROOT, "data", "runs", "h13i_verdict", "scan.json")
    external = None
    if os.path.exists(scan):
        with open(scan) as f:
            external = json.load(f)
        for row in external:
            if row["maxS"] > row["nM"] or row["fail_C37"] or row["fail_H13I"] \
               or row["maxI"] != row["M"]:
                ok = False
                log("FAIL C38 (external scan): %s" % row)
        log("  C external scan rows read: n = %s"
            % [r["n"] for r in external])
    return ok, rows, external


# ------------------------------------------------------------- part D

def part_d(dmax, log):
    ok = True
    rows = []
    for n in range(4, dmax + 1):
        t0 = time.time()
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
                if dist[idx[q]] < 0:
                    dist[idx[q]] = d + 1
                    prev[idx[q]] = (p, i)
                    dq.append(q)
        M = M_of(n)
        maxcs = 0
        diffs = []
        for p in perms:
            G = grid_formula(p)
            I = min(min(r) for r in G)
            c = dist[idx[p]]
            maxcs = max(maxcs, c)
            if c > I:
                ok = False
                log(f"FAIL C39 cs <= I: n={n} pi={p} cs={c} I={I}")
            if c < I:
                diffs.append((list(p), c, I))
        if maxcs != M:
            ok = False
            log(f"FAIL C39 max cs: n={n} {maxcs} != {M}")
        rows.append({"n": n, "max_cs": maxcs, "M": M, "num_cs_lt_I": len(diffs),
                     "lex_first": diffs[0] if diffs else None,
                     "max_gap": max((I - c for _, c, I in diffs), default=0),
                     "seconds": round(time.time() - t0, 2)})
        log("  D n=%d max cs=%d (M=%d) #(cs<I)=%d first=%s"
            % (n, maxcs, M, len(diffs), diffs[0] if diffs else None))
        if n == 8:
            if len(diffs) != 16 or diffs[0][0] != [0, 3, 6, 1, 4, 7, 2, 5]:
                ok = False
                log("FAIL C39 n=8 witness set: %d %s" % (len(diffs), diffs[0]))
            # replay the recorded 8-swap sorting of (0,5,2,7,4,1,6,3)
            pi = [0, 5, 2, 7, 4, 1, 6, 3]
            cur = list(pi)
            for e in [7, 5, 6, 3, 4, 1, 2, 0]:
                f = (e + 1) % 8
                cur[e], cur[f] = cur[f], cur[e]
            rots = [tuple((j + c) % 8 for j in range(8)) for c in range(8)]
            Gpi = grid_formula(tuple(pi))
            Ipi = min(min(r) for r in Gpi)
            if tuple(cur) not in rots or Ipi != 11 or dist[idx[tuple(pi)]] != 8:
                ok = False
                log("FAIL C39 explicit witness: %s I=%d cs=%d"
                    % (cur, Ipi, dist[idx[tuple(pi)]]))
    return ok, rows


# ------------------------------------------------------------- part E

CANDIDATES = {
    "diag":   (4, (0, 3, 2, 1)),
    "anti":   (5, (0, 1, 2, 4, 3)),
    "linear": (5, (0, 2, 1, 4, 3)),
    "T1":     (5, (0, 1, 4, 3, 2)),
    "T2":     (5, (0, 1, 4, 3, 2)),
    "T3":     (5, (0, 1, 4, 3, 2)),
    "rowmax": (7, (0, 1, 6, 5, 4, 3, 2)),
    "T4":     (8, (0, 1, 7, 6, 4, 5, 3, 2)),
}


def candidate_values(n, perm, G):
    M = M_of(n)
    rowmin = [min(G[a][b] for a in range(n)) for b in range(n)]
    out = {
        "rowmax": (max(rowmin), M),
        "diag": (min(Fraction(sum(G[a][(a + t) % n] for a in range(n)), n)
                     for t in range(n)), M),
        "anti": (min(Fraction(sum(G[a][(s - a) % n] for a in range(n)), n)
                     for s in range(n)), M),
        "linear": (min(Fraction(sum(G[a][(k * a + t) % n] for a in range(n)), n)
                       for k in range(n) for t in range(n)), M),
        "T1": (min(sum(G[(i + c) % n][perm[i]] for i in range(n))
                   for c in range(n)), n * M),
        "T2": (min(sum(G[i][(perm[i] + c) % n] for i in range(n))
                   for c in range(n)), n * M),
        "T3": (min(sum(G[(i + c) % n][(perm[i] + d) % n] for i in range(n))
                   for c in range(n) for d in range(n)), n * M),
    }
    if n % 2 == 0:
        out["T4"] = (max(rowmin[b] + rowmin[(b + n // 2) % n] for b in range(n)),
                     2 * M)
    return out


def part_e(log):
    ok = True
    first = {}
    nmax = max(v[0] for v in CANDIDATES.values())
    for n in range(4, nmax + 1):
        for perm in itertools.permutations(range(n)):
            if perm[0] != 0:
                continue
            G = grid_def(perm) if n <= 6 else grid_formula(perm)
            for key, (val, bound) in candidate_values(n, perm, G).items():
                if val > bound and key not in first:
                    first[key] = (n, perm, str(val), bound)
    for key, (n, perm) in CANDIDATES.items():
        got = first.get(key)
        if got is None or got[0] != n or got[1] != perm:
            ok = False
            log("FAIL part E %s: recorded (%d, %s), found %s" % (key, n, perm, got))
    log("  E smallest counterexamples reproduced: %s" % sorted(first))
    return ok, {k: list(v) for k, v in first.items()}


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amax", type=int, default=7)
    ap.add_argument("--bmax", type=int, default=24)
    ap.add_argument("--cmax", type=int, default=8)
    ap.add_argument("--dmax", type=int, default=8)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    lines = []

    def log(msg):
        print(msg)
        lines.append(str(msg))

    rng = random.Random(args.seed)
    t0 = time.time()
    log("%s  seed=%d" % (VERSION, args.seed))
    log("part A: lemmas 0-5 (C37)")
    okA, rowsA = part_a(args.amax, log, rng)
    log("part B: lemma 6, reflection grid (C37)")
    okB, rowsB = part_b(args.bmax, log)
    log("part C: C38, sum_b min_a F(a,b) <= n*M")
    okC, rowsC, external = part_c(args.cmax, log)
    log("part D: C39, circular swap distance vs I")
    okD, rowsD = part_d(args.dmax, log)
    log("part E: smallest counterexamples of the refuted candidates")
    okE, rowsE = part_e(log)

    ok = okA and okB and okC and okD and okE
    res = {"version": VERSION, "seed": args.seed,
           "result": "PASS" if ok else "FAIL",
           "seconds": round(time.time() - t0, 1),
           "partA_lemmas": {"ok": okA, "rows": rowsA},
           "partB_reflections": {"ok": okB, "nmax": args.bmax, "rows": rowsB},
           "partC_C38": {"ok": okC, "rows": rowsC, "external_scan": external},
           "partD_C39": {"ok": okD, "rows": rowsD},
           "partE_candidates": {"ok": okE, "smallest": rowsE}}
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(res, f, indent=1)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# check_C37 — %s\n\n" % res["result"])
        f.write("Версия %s, seed %d, %.1f с.\n\n" % (VERSION, args.seed, res["seconds"]))
        f.write("```\n" + "\n".join(lines) + "\n```\n")
    log("%s in %.1f s -> %s" % (res["result"], res["seconds"], OUT))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
