"""Toric class of a permutation: exact identity inv = D2/n + const, mean of inv,
and averaged variants of H13-I (session 9).

Model (docs/notes/h13_line_model.md sec. 0).  For pi: Z_n -> Z_n and a double cut
(a, b) in Z_n^2 the line is w_j = (pi(j + a) - b) mod n, j = 0..n-1; equivalently
the point i of the torus gets coordinates X_i = (i - a) mod n, Y_i = (pi(i) - b) mod n.
  inv(a, b) = number of inversions of w;
  D2(a, b)  = sum_j (j - w_j)^2 = sum_i (X_i - Y_i)^2   (squared Spearman);
  I(pi)     = min over the n^2 cuts of inv(a, b)        (H13-I: I <= floor((n-1)^2/4));
  A(pi)     = mean over the n^2 cuts of inv(a, b).

Measured here (exhaustively over all pi unless stated):
  (1) inv(a, b) - D2(a, b)/n is the same for all n^2 cuts (theorem A, C37) and equals
      (n^2-1)/12 - W/n^2 with W = sum_{i,j} f(j-i) f(pi(j)-pi(i)), f(t) = t - (n-1)/2;
  (2) max_pi A(pi) = (n-1)(2n-1)/6, attained exactly on the n reflections (theorem B,
      C38); hence I(pi) <= floor((n-1)(2n-1)/6) for every pi;
  (3) mean over a of min over b of inv(a, b) <= floor((n-1)^2/4), equality exactly on
      the reflections (observation C39, finite data only);
  (4) counterexamples to the averaged variants over "lines" of cuts
      {(a, mu*a + c)} and over columns (rejected family);
  (5) the second moment refinement A - 6 Var/(n^2-1) of (2) on reflections: how far
      it stays above floor((n-1)^2/4).
Usage: python3 experiments/toric_identity.py --nmax 8 [--perms 9]
Output: data/runs/toric_identity/raw*.md, report*.json.  Version toric_identity-1.0.
"""

import argparse
import itertools
import json
import os
import platform
import sys
import time
from fractions import Fraction

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
VERSION = "toric_identity-1.0"
OUT = os.path.join(ROOT, "data", "runs", "toric_identity")


def M(n):
    return ((n - 1) ** 2) // 4


def inv_brute(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def inv_table(pi):
    """inv(a, b) for all n^2 cuts, by the two step recurrences (O(n^2))."""
    n = len(pi)
    ipi = [0] * n
    for i, v in enumerate(pi):
        ipi[v] = i
    tab = [[0] * n for _ in range(n)]
    cur = inv_brute(list(pi))
    tab[0][0] = cur
    for b in range(n - 1):
        cur += (n - 1) - 2 * ipi[b]          # value cut: min goes to max
        tab[0][b + 1] = cur
    for b in range(n):
        cur = tab[0][b]
        row = tab
        for a in range(n - 1):
            cur += (n - 1) - 2 * ((pi[a] - b) % n)   # position cut: first goes last
            row[a + 1][b] = cur
    return tab


def d2_table(pi):
    n = len(pi)
    tab = [[0] * n for _ in range(n)]
    for a in range(n):
        for b in range(n):
            tab[a][b] = sum((j - (pi[(j + a) % n] - b) % n) ** 2 for j in range(n))
    return tab


def const_sawtooth(pi):
    """(n^2-1)/12 - W/n^2 in exact arithmetic."""
    n = len(pi)
    f = [Fraction(2 * t - (n - 1), 2) for t in range(n)]
    W = sum(f[(j - i) % n] * f[(pi[j] - pi[i]) % n] for i in range(n) for j in range(n))
    return Fraction(n * n - 1, 12) - Fraction(W, n * n)


def is_reflection(pi):
    n = len(pi)
    s = (pi[0] + 0) % n
    return all((pi[i] + i) % n == s for i in range(n))


def run_n(n, identity_sample, skip_lines=False):
    """Exhaustive pass over S_n."""
    t0 = time.time()
    res = {"n": n, "M": M(n)}
    sumA_max = -1          # max over pi of n^2 * A(pi)
    argA = []
    max_avgmin = -1
    arg_avgmin = []
    max_I = -1
    bad_line = 0
    worst_line = None
    checked_identity = 0
    npi = 0
    for pi in itertools.permutations(range(n)):
        npi += 1
        tab = inv_table(pi)
        # (1) identity, on a sample of permutations (D2 table costs O(n^3))
        if npi <= identity_sample or npi % 1009 == 0:
            checked_identity += 1
            dt = d2_table(pi)
            vals = {Fraction(tab[a][b]) - Fraction(dt[a][b], n)
                    for a in range(n) for b in range(n)}
            assert len(vals) == 1, ("identity broken", pi, sorted(vals)[:4])
            assert vals.pop() == const_sawtooth(pi), ("const formula", pi)
        # (2) mean
        s = sum(sum(r) for r in tab)
        if s > sumA_max:
            sumA_max, argA = s, [pi]
        elif s == sumA_max:
            argA.append(pi)
        I = min(min(r) for r in tab)
        max_I = max(max_I, I)
        # (3) mean over a of min over b
        am = sum(min(r) for r in tab)          # n * (mean over a of min over b)
        if am > max_avgmin:
            max_avgmin, arg_avgmin = am, [pi]
        elif am == max_avgmin:
            arg_avgmin.append(pi)
        # (4) line families
        if skip_lines:
            continue
        best = None                            # n * (mean of inv along the family)
        for mu in range(n):
            for c in range(n):
                v = sum(tab[a][(mu * a + c) % n] for a in range(n))
                if best is None or v < best:
                    best = v
        for a in range(n):
            v = sum(tab[a])
            if v < best:
                best = v
        if best > n * M(n):
            bad_line += 1
            if worst_line is None or best > worst_line[0]:
                worst_line = (best, pi)
    res["perms"] = npi
    res["identity_checked"] = checked_identity
    res["max_n2A"] = sumA_max
    res["A_bound_n2"] = n * n * (n - 1) * (2 * n - 1) // 6
    res["argmax_A_count"] = len(argA)
    res["argmax_A_all_reflections"] = all(is_reflection(p) for p in argA)
    res["max_I"] = max_I
    res["I_general_bound"] = (n - 1) * (2 * n - 1) // 6
    res["max_avg_a_min_b"] = str(Fraction(max_avgmin, n))
    res["avgmin_tight_count"] = len(arg_avgmin)
    res["avgmin_tight_all_reflections"] = all(is_reflection(p) for p in arg_avgmin)
    res["line_family_failures"] = "не считалось" if skip_lines else bad_line
    res["line_family_worst"] = ((str(Fraction(worst_line[0], n)), list(worst_line[1]))
                                if worst_line else None)
    res["seconds"] = round(time.time() - t0, 1)
    return res


def second_moment_on_reflections(nmax):
    """A - 6 Var/(n^2-1) on the reflection pi(i) = -i, against floor((n-1)^2/4)."""
    out = []
    for n in range(4, nmax + 1):
        pi = tuple((-i) % n for i in range(n))
        tab = inv_table(pi)
        vals = [tab[a][b] for a in range(n) for b in range(n)]
        A = Fraction(sum(vals), n * n)
        var = Fraction(sum(v * v for v in vals), n * n) - A * A
        bound = A - 6 * var / (n * n - 1)
        out.append({"n": n, "A": str(A), "var": str(var), "refined": str(bound),
                    "float": round(float(bound), 3), "M": M(n), "I": min(vals)})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--nmin", type=int, default=4)
    ap.add_argument("--identity_sample", type=int, default=50)
    ap.add_argument("--skip_lines", action="store_true")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    rows = [run_n(n, args.identity_sample, args.skip_lines)
            for n in range(args.nmin, args.nmax + 1)]
    refine = second_moment_on_reflections(min(args.nmax + 4, 20))
    report = {
        "version": VERSION,
        "python": platform.python_version(),
        "args": vars(args),
        "rows": rows,
        "second_moment_refinement_on_reflections": refine,
        "seconds": round(time.time() - t0, 1),
    }
    suffix = "" if args.nmin == 4 else f"_n{args.nmin}"
    with open(os.path.join(OUT, f"report{suffix}.json"), "w") as fh:
        json.dump(report, fh, indent=1)
    lines = [f"# toric_identity ({VERSION}) n = {args.nmin}..{args.nmax}", ""]
    for r in rows:
        ok_mean = r["max_n2A"] == r["A_bound_n2"]
        lines += [
            f"n = {r['n']} ({r['perms']} перестановок, {r['seconds']} с)",
            f"  тождество inv - D2/n = const и формула const: проверено на "
            f"{r['identity_checked']} pi (все n^2 разрезов, инверсии перебором)",
            f"  max n^2 A(pi) = {r['max_n2A']}, граница n^2 (n-1)(2n-1)/6 = "
            f"{r['A_bound_n2']} -> {'равенство' if ok_mean else 'НАРУШЕНО'}; "
            f"argmax: {r['argmax_A_count']} шт., все отражения: {r['argmax_A_all_reflections']}",
            f"  max I = {r['max_I']} <= floor((n-1)(2n-1)/6) = {r['I_general_bound']} "
            f"(и <= M_n = {r['M']})",
            f"  max avg_a min_b inv = {r['max_avg_a_min_b']} (M_n = {r['M']}), "
            f"достигается на {r['avgmin_tight_count']} pi, все отражения: "
            f"{r['avgmin_tight_all_reflections']}",
            f"  семейства-прямые: минимум среднего > M_n у {r['line_family_failures']} pi, "
            f"худший {r['line_family_worst']}",
            "",
        ]
    lines += ["Уточнение вторым моментом на отражениях (A - 6Var/(n^2-1) против M_n):", ""]
    for r in refine:
        lines.append(f"  n = {r['n']:3d}: A - 6Var/(n^2-1) = {r['float']:9.3f}, "
                     f"M_n = {r['M']:5d}, I = {r['I']:5d}")
    lines += ["", f"Всего {report['seconds']} с."]
    with open(os.path.join(OUT, f"raw{suffix}.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
