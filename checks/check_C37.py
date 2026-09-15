"""Checker for C37 (session 9): exactness of the constant in H13-I and the
obstructions to averaging over families of double cuts.

Independent implementation: every line is built explicitly and its inversions
are counted by enumerating pairs; the incremental table of experiments/
toric_cuts.py is NOT used here.

  Part A (Theorem 1, docs/proofs/C37_toric_reflection.md):
    for every reflection pi_h(i) = (h - i) mod n and EVERY double cut (p, s)
    inv(line) = C(n,2) - (t+1)(n-1-t) with t = (h - p - s) mod n, hence
    inv >= floor((n-1)^2/4), and I(pi_h) = floor((n-1)^2/4).
    All h, all n^2 cuts, 4 <= n <= AMAX_A (default 20).
  Part B (Lemmas 3-5): the triple statistic T(pi) is invariant under the two
    generators of the toric action (hence under all n^2 shifts), (n-2)*I >= T,
    and T = C(n,3) exactly on the reflections -- exhaustively for 4 <= n <= 7.
    Finite observation: {pi : I(pi) = floor((n-1)^2/4)} is exactly the set of
    the n reflections, exhaustively for 4 <= n <= 8.
  Part C (Lemma 6): the per-pair cut count n*dx + n*dy - 2*dx*dy and the mean
    formula mean inv = C(n,2) - (1/n^2) sum_{ordered pairs} dx*dy
    (exhaustive 4 <= n <= 6, sample at n = 7, 8); the reflection mean equals
    C(n,2) - (n^2-1)/6 and exceeds floor((n-1)^2/4) by at least (n^2-1)/12
    (exact arithmetic, 4 <= n <= 200).
  Part D (affine no-go): for gcd(a, n) = 1 and pi(i) = a*i + b, inv(p, s)
    depends only on r = (a*p + b - s) mod n, each r comes from exactly n cuts,
    and the point-anchored family F(alpha, beta) attains I(pi)
    (4 <= n <= 12); outside the affine inputs that family fails
    (witnesses n = 5, 7, 8).

Usage: python3 checks/check_C37.py [--amax-a 20] [--amax-eq 8]
Output: data/runs/check_C37/report.json, report.md.
"""

import argparse
import itertools
import json
import os
import random
import sys
import time
from fractions import Fraction

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
VERSION = "check_C37-1.0"
OUT = os.path.join(ROOT, "data", "runs", "check_C37")


def bound(n):
    return ((n - 1) ** 2) // 4


def C2(n):
    return n * (n - 1) // 2


def line(pi, n, p, s):
    return [(pi[(p + j) % n] - s) % n for j in range(n)]


def inv_of(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def I_direct(pi, n, stop_below=None):
    """min over all n^2 cuts of inv; with stop_below returns as soon as a cut
    with inv < stop_below is found (returning that value)."""
    best = None
    for p in range(n):
        for s in range(n):
            v = inv_of(line(pi, n, p, s))
            if best is None or v < best:
                best = v
                if stop_below is not None and best < stop_below:
                    return best
    return best


def triples_T(pi, n):
    T = 0
    for a, b, c in itertools.combinations(range(n), 3):
        x, y, z = pi[a], pi[b], pi[c]
        k = (x > y) + (x > z) + (y > z)
        if k % 2 == 1:
            T += 1
    return T


def reflections(n):
    return [tuple((h - i) % n for i in range(n)) for h in range(n)]


def is_reflection(pi, n):
    return all((pi[(i + 1) % n] - pi[i]) % n == n - 1 for i in range(n))


def affine_perms(n):
    out = []
    for a in range(1, n):
        if all((a * i) % n for i in range(1, n)):
            for b in range(n):
                out.append((a, b, tuple((a * i + b) % n for i in range(n))))
    return out


def part_a(amax, log):
    ok = True
    for n in range(4, amax + 1):
        t0 = time.time()
        M = bound(n)
        for h in range(n):
            pi = tuple((h - i) % n for i in range(n))
            mn = None
            for p in range(n):
                for s in range(n):
                    t = (h - p - s) % n
                    v = inv_of(line(pi, n, p, s))
                    if v != C2(n) - (t + 1) * (n - 1 - t):
                        ok = False
                        log(f"FAIL part A formula: n={n} h={h} p={p} s={s}")
                    if v < M:
                        ok = False
                        log(f"FAIL part A lower bound: n={n} h={h} p={p} s={s} inv={v} < {M}")
                    mn = v if mn is None else min(mn, v)
            if mn != M:
                ok = False
                log(f"FAIL part A minimum: n={n} h={h} I={mn} != {M}")
        log(f"part A n={n}: all {n} reflections, all {n*n} cuts: "
            f"inv = C(n,2)-(t+1)(n-1-t) >= {M}, I = {M}, {time.time()-t0:.1f} s")
    return ok


def part_b(amax_t, amax_eq, log):
    ok = True
    for n in range(4, amax_t + 1):
        t0 = time.time()
        cn3 = n * (n - 1) * (n - 2) // 6
        for pi in itertools.permutations(range(n)):
            T = triples_T(pi, n)
            sh_p = tuple(pi[(i + 1) % n] for i in range(n))
            sh_s = tuple((pi[i] + 1) % n for i in range(n))
            if triples_T(sh_p, n) != T or triples_T(sh_s, n) != T:
                ok = False
                log(f"FAIL part B invariance: n={n} pi={pi}")
            I = I_direct(pi, n)
            if (n - 2) * I < T:
                ok = False
                log(f"FAIL part B (n-2)I>=T: n={n} pi={pi} I={I} T={T}")
            if (T == cn3) != is_reflection(pi, n):
                ok = False
                log(f"FAIL part B T=C(n,3) <=> reflection: n={n} pi={pi} T={T}")
        log(f"part B n={n}: exhaustive, T invariant under both generators, "
            f"(n-2)I >= T, T = C(n,3) exactly on reflections, {time.time()-t0:.1f} s")
    for n in range(4, amax_eq + 1):
        t0 = time.time()
        M = bound(n)
        eq = []
        for pi in itertools.permutations(range(n)):
            v = I_direct(pi, n, stop_below=M)
            if v > M:
                ok = False
                log(f"FAIL part B H13-I violated: n={n} pi={pi} I={v}")
            elif v == M:
                eq.append(pi)
        if sorted(eq) != sorted(reflections(n)):
            ok = False
            log(f"FAIL part B equality set: n={n} size={len(eq)}")
        log(f"part B n={n}: exhaustive, I <= {M} and equality exactly on the "
            f"{len(eq)} reflections, {time.time()-t0:.1f} s")
    return ok


def part_c(log, rng):
    ok = True
    for n in range(4, 9):
        if n <= 6:
            perms = list(itertools.permutations(range(n)))
            cover = "exhaustive"
        else:
            perms = [tuple(rng.sample(range(n), n)) for _ in range(200)] + reflections(n)
            cover = f"sample of {len(perms)}"
        for pi in perms:
            pts = [(i, pi[i]) for i in range(n)]
            tot = 0
            for (x1, y1), (x2, y2) in itertools.combinations(pts, 2):
                dx = (x2 - x1) % n
                dy = (y2 - y1) % n
                cnt = sum(1 for p in range(n) for s in range(n)
                          if (((x1 - p) % n < (x2 - p) % n) !=
                              ((y1 - s) % n < (y2 - s) % n)))
                if cnt != n * dx + n * dy - 2 * dx * dy:
                    ok = False
                    log(f"FAIL part C pair count: n={n} pi={pi}")
                tot += cnt
            real = sum(inv_of(line(pi, n, p, s)) for p in range(n) for s in range(n))
            if real != tot:
                ok = False
                log(f"FAIL part C mean: n={n} pi={pi} {real} != {tot}")
            sdxdy = 0
            for (x1, y1), (x2, y2) in itertools.permutations(pts, 2):
                sdxdy += ((x2 - x1) % n) * ((y2 - y1) % n)
            if Fraction(real, n * n) != C2(n) - Fraction(sdxdy, n * n):
                ok = False
                log(f"FAIL part C mean formula: n={n} pi={pi}")
        log(f"part C n={n}: {cover}: per-pair cut count and mean formula hold")
    for n in range(4, 201):
        mean = Fraction(C2(n)) - Fraction(n * n - 1, 6)
        pi = tuple((-i) % n for i in range(n))
        if n <= 9:
            real = Fraction(sum(inv_of(line(pi, n, p, s))
                                for p in range(n) for s in range(n)), n * n)
            if real != mean:
                ok = False
                log(f"FAIL part C reflection mean: n={n}")
        if mean - bound(n) < Fraction(n * n - 1, 12):
            ok = False
            log(f"FAIL part C gap: n={n}")
    log("part C: reflection mean = C(n,2) - (n^2-1)/6 (4 <= n <= 9 direct), "
        "mean - floor((n-1)^2/4) >= (n^2-1)/12 in exact arithmetic (4 <= n <= 200)")
    return ok


def part_d(amax, log):
    ok = True
    for n in range(4, amax + 1):
        aff = affine_perms(n)
        for a, b, pi in aff:
            byr = {}
            for p in range(n):
                for s in range(n):
                    r = (a * p + b - s) % n
                    byr.setdefault(r, []).append(inv_of(line(pi, n, p, s)))
            if len(byr) != n or any(len(v) != n or len(set(v)) != 1 for v in byr.values()):
                ok = False
                log(f"FAIL part D affine: n={n} a={a} b={b}")
            # point-anchored family attains I
            I = min(v[0] for v in byr.values())
            best = None
            for al in range(n):
                for be in range(n):
                    tot = sum(inv_of(line(pi, n, (i + al) % n, (pi[i] + be) % n))
                              for i in range(n))
                    best = tot if best is None else min(best, tot)
            if best != n * I:
                ok = False
                log(f"FAIL part D anchored on affine: n={n} a={a} b={b} {best} vs {n*I}")
        log(f"part D n={n}: {len(aff)} affine inputs: inv depends only on r, "
            f"each r from exactly n cuts, anchored family attains I(pi)")
    for n, pi in ((5, (0, 1, 4, 3, 2)), (7, (0, 2, 1, 6, 5, 4, 3)),
                  (8, (5, 7, 6, 4, 3, 0, 2, 1))):
        best = None
        for al in range(n):
            for be in range(n):
                tot = sum(inv_of(line(pi, n, (i + al) % n, (pi[i] + be) % n))
                          for i in range(n))
                best = tot if best is None else min(best, tot)
        I = I_direct(pi, n)
        good = best > n * bound(n)
        if not good:
            ok = False
        log(f"part D witness n={n} pi={pi}: best anchored family average = "
            f"{Fraction(best, n)} > floor((n-1)^2/4) = {bound(n)} -> "
            f"{'family insufficient (as claimed)' if good else 'FAIL'} "
            f"(true I(pi) = {I})")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amax-a", type=int, default=20)
    ap.add_argument("--amax-t", type=int, default=7)
    ap.add_argument("--amax-eq", type=int, default=8)
    ap.add_argument("--amax-d", type=int, default=12)
    ap.add_argument("--seed", type=int, default=20260915)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    rng = random.Random(args.seed)
    lines = []

    def log(msg):
        print(msg)
        lines.append(msg)

    log(f"== {VERSION} args={vars(args)}")
    t0 = time.time()
    ok = part_a(args.amax_a, log)
    ok = part_b(args.amax_t, args.amax_eq, log) and ok
    ok = part_c(log, rng) and ok
    ok = part_d(args.amax_d, log) and ok
    log(f"verdict: {'PASS' if ok else 'FAIL'} ({time.time()-t0:.1f} s)")
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump({"version": VERSION, "args": vars(args),
                   "verdict": "PASS" if ok else "FAIL", "log": lines}, f, indent=1)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# check_C37 — точность константы H13-I и запреты усреднений\n\n")
        f.write(f"Команда: `python3 checks/check_C37.py --amax-a {args.amax_a} "
                f"--amax-t {args.amax_t} --amax-eq {args.amax_eq} "
                f"--amax-d {args.amax_d}`. Версия: {VERSION}.\n\n```text\n")
        f.write("\n".join(lines))
        f.write("\n```\n")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
