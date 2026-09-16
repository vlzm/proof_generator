"""Checker for C37 (session 9): concordance-sum identity and averaging insufficiency.

Double cut (q, c): position origin a = q + 1, value origin b = q + 1 - c (both
mod n); w_j = pi(a + j) - b mod n, j = 0..n-1 (matches experiments/line_model.py
and experiments/line_profile.c).  (q, c) <-> (a, b) is a bijection of Z_n x Z_n,
so I(pi) = min over (q, c) of inv(w) = min over (a, b) of inv(w).

Part A (identity): for a fixed pi, sum over all n^2 cuts (a, b) of noninv(w)
equals sum over unordered pairs {i1, i2} (i1 < i2 as plain integers) of
  n^2 - n*(di + dj) + 2*di*dj,
where di = i2 - i1 (in 1..n-1) and dj = (pi(i2) - pi(i1)) mod n (in 1..n-1).
Verified against brute force over all n^2 cuts, all permutations, 2 <= n <= AMAX.

Part B (averaging insufficiency, exact): for sigma_n (pi(i) = (1 - i) mod n)
and every rotation L^k sigma_n, dj = n - di identically for every pair, so the
per-cut average of noninv(w) over the n^2 cuts equals exactly (n^2 - 1)/6,
a Theta(n^2) shortfall below the H13-I target floor(n^2/4) (not O(n)): no
per-pair-independent (weighted or uniform) averaging bound over the n^2 cuts
can establish H13-I, even though sigma_n itself meets the target exactly at
one specific cut (C32/C33).  Checked in exact rational arithmetic, no
permutation enumeration needed, for 4 <= n <= NMAX.

Usage: python3 checks/check_C37.py [--amax 7] [--nmax 200]
Output: data/runs/check_C37/report.json, report.md.
"""

import argparse
import itertools
import json
import os
import time
from fractions import Fraction

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

VERSION = "check_C37-1.0"
OUT = os.path.join(ROOT, "data", "runs", "check_C37")


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def noninv_bruteforce_sum(pi):
    """Sum over all n^2 cuts (a, b) of noninv(w); w_j = pi[(a+j) % n] - b mod n."""
    n = len(pi)
    total = 0
    for a in range(n):
        for b in range(n):
            w = [(pi[(a + j) % n] - b) % n for j in range(n)]
            total += n * (n - 1) // 2 - inversions(w)
    return total


def noninv_formula_sum(pi):
    n = len(pi)
    s = 0
    for i1 in range(n):
        for i2 in range(i1 + 1, n):
            di = i2 - i1
            dj = (pi[i2] - pi[i1]) % n
            s += n * n - n * (di + dj) + 2 * di * dj
    return s


def part_a(amax, log):
    ok = True
    checked = 0
    for n in range(2, amax + 1):
        t0 = time.time()
        for pi in itertools.permutations(range(n)):
            bf = noninv_bruteforce_sum(list(pi))
            fm = noninv_formula_sum(list(pi))
            checked += 1
            if bf != fm:
                ok = False
                log(f"FAIL part A: n={n} pi={pi} bruteforce={bf} formula={fm}")
        log(f"part A: n={n} all permutations OK ({time.time()-t0:.1f}s)")
    return ok, checked


def part_b(nmax, log):
    """Exact average of noninv over the n^2 cuts for sigma_n (and its rotations),
    computed from the closed form (dj = n - di for every pair): should equal
    (n^2 - 1) / 6 exactly, and fall short of target floor(n^2/4) by Theta(n^2)."""
    ok = True
    rows = []
    for n in range(4, nmax + 1):
        sigma = [(1 - i) % n for i in range(n)]
        s = noninv_formula_sum(sigma)
        avg = Fraction(s, n * n)
        expected = Fraction(n * n - 1, 6)
        target = (n - 1) * (n - 1) // 4
        gap = float(avg - target)
        row = {
            "n": n,
            "avg": float(avg),
            "expected_(n2-1)/6": float(expected),
            "target_floor((n-1)^2/4)_noninv_equiv": (n * n) // 4,
            "gap_avg_minus_target": gap,
        }
        rows.append(row)
        if avg != expected:
            ok = False
            log(f"FAIL part B: n={n} avg={avg} expected={expected}")
        # also check a rotation L^k sigma_n has the same average (dj = n - di still)
        for k in (1, n // 3, n // 2):
            rot = [sigma[(i + k) % n] for i in range(n)]
            s2 = noninv_formula_sum(rot)
            if Fraction(s2, n * n) != expected:
                ok = False
                log(f"FAIL part B rotation: n={n} k={k}")
    return ok, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amax", type=int, default=7)
    ap.add_argument("--nmax", type=int, default=200)
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    lines = []

    def log(msg):
        print(msg)
        lines.append(msg)

    log(f"{VERSION} amax={args.amax} nmax={args.nmax}")
    okA, checkedA = part_a(args.amax, log)
    okB, rowsB = part_b(args.nmax, log)
    ok = okA and okB
    log(f"RESULT: {'PASS' if ok else 'FAIL'}")

    report = {
        "version": VERSION,
        "amax": args.amax,
        "nmax": args.nmax,
        "part_a_pass": okA,
        "part_a_permutations_checked": checkedA,
        "part_b_pass": okB,
        "part_b_rows": rowsB,
        "overall_pass": ok,
    }
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(report, f, indent=2)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write(f"# {VERSION} report\n\n")
        f.write(f"Part A (identity, brute force vs formula): {'PASS' if okA else 'FAIL'}, "
                f"{checkedA} permutations, 2 <= n <= {args.amax}.\n\n")
        f.write(f"Part B (sigma_n average = (n^2-1)/6, exact): {'PASS' if okB else 'FAIL'}, "
                f"4 <= n <= {args.nmax}.\n\n")
        f.write("| n | avg | (n^2-1)/6 | target floor(n^2/4) | gap |\n|---|---|---|---|---|\n")
        for r in rowsB:
            f.write(f"| {r['n']} | {r['avg']:.4f} | {r['expected_(n2-1)/6']:.4f} | "
                    f"{r['target_floor((n-1)^2/4)_noninv_equiv']} | {r['gap_avg_minus_target']:.4f} |\n")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
