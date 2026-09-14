"""h13i_leads.py -- reproducible record of the approaches to H13-I tried in session 9.

H13-I:  I(pi) = min over the n^2 double cuts (a,b) of inv(w), w_j = (pi(a+j)-b) mod n,
satisfies I(pi) <= floor((n-1)^2/4) for every toroidal permutation pi of Z_n.

Every lead below is either REFUTED (a counterexample is printed) or reported as
finite data.  Notation follows docs/proofs/C37_inversion_spearman.md:

    psi(t) = (t mod n) - (n-1)/2,   V(a,b) = sum_i psi(i-a) psi(pi(i)-b),
    inv(a,b) = inv_avg(pi) - (2/n) V(a,b),      Q(pi) = sum_i V(i,pi(i)),
    c_n = (n^2-1)/4 - floor((n-1)^2/4),
    H13-I(pi) <=> 2 V*(pi) + Q(pi)/n >= n c_n          (C37.7)
    3 V*(pi) >= n c_n  =>  H13-I(pi)                   (C37.8)

Leads:
  L0  averaging over all n^2 cuts                          (refuted, session 8)
  L1  averaging over a diagonal {(a, a+s)} / antidiagonal {(a, s-a)}
  L2  averaging over the point-cut family {(i, pi(i)+r)} and {(i+q, pi(i)+r)}
  L3  Mantel route (a): a cut whose line is 321-avoiding (triangle-free
      inversion graph) would give inv <= floor(n^2/4) by Mantel -- too weak by
      floor(n/2), and such a cut need not exist
  L4  Mantel route (b): a cut with a balanced block split (the first ceil(n/2)
      or floor(n/2) line entries are exactly the smallest values) would give the
      exact bound -- such a cut need not exist
  L5  sufficient condition (C37.8): does V* >= n c_n / 3 hold for every pi?
  L6  two explicit lower bounds on V*:
        L6a (diagonal shifts)  V* >= (1/n) max_r sum_i C(pi(i)-i-r), C = autocorr(psi)
        L6b (point cuts)       V* >= (1/n) max_r sum_i V(i, pi(i)+r)

Usage: python3 experiments/h13i_leads.py [--nmax 8]
Output: data/runs/h13i_leads/report.json, report.md.  Version h13i_leads-1.0.
"""

import argparse
import itertools
import json
import os
import sys
import time
from fractions import Fraction as Fr

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "data", "runs", "h13i_leads")
VERSION = "h13i_leads-1.0"


def inv_table(pi):
    """T[a][b] = number of inversions of the line at cut (a,b) (O(n^2) total)."""
    n = len(pi)
    T = [[0] * n for _ in range(n)]
    T[0][0] = sum(1 for i in range(n) for j in range(i + 1, n) if pi[i] > pi[j])
    posval = [0] * n
    for i, v in enumerate(pi):
        posval[v] = i
    for a in range(n):
        if a > 0:
            T[a][0] = T[a - 1][0] + (n - 1) - 2 * pi[a - 1]
        for b in range(1, n):
            T[a][b] = T[a][b - 1] + (n - 1) - 2 * ((posval[b - 1] - a) % n)
    return T


def V4_table(pi):
    """4*V(a,b)."""
    n = len(pi)
    return [[sum((2 * ((i - a) % n) - (n - 1)) * (2 * ((pi[i] - b) % n) - (n - 1))
                 for i in range(n)) for b in range(n)] for a in range(n)]


def reps(n):
    for tail in itertools.permutations(range(1, n)):
        yield (0,) + tail


def has_321(w):
    n = len(w)
    for i in range(n):
        for j in range(i + 1, n):
            if w[j] < w[i]:
                for k in range(j + 1, n):
                    if w[k] < w[j]:
                        return True
    return False


def run_n(n):
    bound = ((n - 1) ** 2) // 4
    c4n2 = n * n * (n * n - 1) - 4 * n * n * bound       # = 4 n^2 c_n
    res = {"n": n, "bound": bound, "c_n": str(Fr(n * n - 1, 4) - bound)}
    fails = {k: [] for k in ("L1diag", "L1anti", "L2r", "L2qr", "L3", "L4", "L5", "L6")}
    cnt = 0
    ks = {n // 2, (n + 1) // 2}
    for pi in reps(n):
        cnt += 1
        T = inv_table(pi)
        Vt = V4_table(pi)
        I = min(min(row) for row in T)
        # ---- L1: one-parameter averaging over (anti)diagonals
        if min(sum(T[a][(a + s) % n] for a in range(n)) for s in range(n)) > n * bound:
            fails["L1diag"].append(pi)
        if min(sum(T[a][(s - a) % n] for a in range(n)) for s in range(n)) > n * bound:
            fails["L1anti"].append(pi)
        # ---- L2: point-cut families
        if min(sum(T[i][(pi[i] + r) % n] for i in range(n)) for r in range(n)) > n * bound:
            fails["L2r"].append(pi)
        if min(sum(T[(i + q) % n][(pi[i] + r) % n] for i in range(n))
               for q in range(n) for r in range(n)) > n * bound:
            fails["L2qr"].append(pi)
        # ---- L3 / L4: Mantel routes
        ok321 = False
        okblk = False
        for a in range(n):
            for b in range(n):
                w = [(pi[(a + j) % n] - b) % n for j in range(n)]
                if not ok321 and not has_321(w):
                    ok321 = True
                if not okblk:
                    for k in ks:
                        if max(w[:k]) < k:
                            okblk = True
                            break
            if ok321 and okblk:
                break
        if not ok321:
            fails["L3"].append(pi)
        if not okblk:
            fails["L4"].append(pi)
        # ---- L5 / L6
        V4s = max(max(row) for row in Vt)
        Q4 = sum(Vt[i][pi[i]] for i in range(n))
        assert (2 * n * V4s + Q4 >= c4n2) == (I <= bound)          # C37.7, re-checked
        if 3 * n * V4s < c4n2:                                     # C37.8 fails
            fails["L5"].append(pi)
        # explicit lower bounds on V* (in 4x integers)
        C4 = [sum((2 * t - (n - 1)) * (2 * ((t + m) % n) - (n - 1)) for t in range(n))
              for m in range(n)]
        L6a = max(Fr(sum(C4[(pi[i] - i - r) % n] for i in range(n)), n) for r in range(n))
        L6b = max(Fr(sum(Vt[i][(pi[i] + r) % n] for i in range(n)), n) for r in range(n))
        if 3 * n * max(L6a, L6b) < c4n2:
            fails["L6"].append(pi)
    res["representatives"] = cnt
    for k, v in fails.items():
        res[k] = {"failures": len(v), "example": list(v[0]) if v else None}
    return res


def families(n, seed=20260914):
    """Structured and random pi for the large-n sampling ladder (AGENTS rule 9)."""
    import random
    from math import gcd
    rng = random.Random(seed + n)
    out = [("id", list(range(n))), ("reflection", [(-i) % n for i in range(n)]),
           ("rev_block", [(n - 1 - i) % n for i in range(n)])]
    for g in range(2, n):
        if gcd(g, n) == 1:
            out.append((f"affine g={g}", [(g * i) % n for i in range(n)]))
    for k in (2, 3, n // 4 or 1):
        if 0 < k < n:
            out.append((f"blocks k={k}", [((i % k) * (n // k) + i // k) % n for i in range(n)]))
    for t in range(5):
        p = list(range(n))
        rng.shuffle(p)
        out.append((f"random {t}", p))
    return out


def run_families(ns):
    rows = []
    for n in ns:
        bound = ((n - 1) ** 2) // 4
        c4n2 = n * n * (n * n - 1) - 4 * n * n * bound
        worst_h13 = None
        worst_suff = None
        for name, pi in families(n):
            T = inv_table(pi)
            I = min(min(r) for r in T)
            Vt = V4_table(pi)
            V4s = max(max(r) for r in Vt)
            # margins, normalised to integers: H13-I margin = bound - I
            if worst_h13 is None or bound - I < worst_h13[0]:
                worst_h13 = (bound - I, name)
            s = 3 * n * V4s - c4n2                     # >= 0 means C37.8 holds
            if worst_suff is None or s < worst_suff[0]:
                worst_suff = (s, name)
        rows.append({"n": n, "bound": bound,
                     "min_margin_H13I": worst_h13[0], "at": worst_h13[1],
                     "min_margin_C37.8_x4n": worst_suff[0], "at2": worst_suff[1]})
        print(f"n={n}: min (bound - I) = {worst_h13[0]} at {worst_h13[1]}; "
              f"min 4n-scaled slack of C37.8 = {worst_suff[0]} at {worst_suff[1]}")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--families", type=str, default="20,50,64,100,101",
                    help="comma separated n for the structured/random sampling ladder")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    rows = [run_n(n) for n in range(4, args.nmax + 1)]
    fam_ns = [int(x) for x in args.families.split(",") if x.strip()]
    fam = run_families(fam_ns) if fam_ns else []
    dt = time.time() - t0
    rep = {"version": VERSION, "seconds": round(dt, 1),
           "coverage": f"all pi with pi(0)=0 (a transversal of the torus classes), 4 <= n <= {args.nmax}"
                       f"; structured/random families at n = {fam_ns}",
           "rows": rows, "families": fam}
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(rep, f, indent=1, ensure_ascii=False)
    hdr = ["n", "bound", "L1diag", "L1anti", "L2r", "L2qr", "L3", "L4", "L5", "L6"]
    lines = ["| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)]
    for r in rows:
        lines.append("| " + " | ".join(
            [str(r["n"]), str(r["bound"])] +
            [str(r[k]["failures"]) for k in hdr[2:]]) + " |")
    ex = []
    for r in rows:
        for k in hdr[2:]:
            if r[k]["example"]:
                ex.append(f"- n={r['n']} {k}: pi = {tuple(r[k]['example'])}")
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write(f"# h13i_leads ({VERSION})\n\nЧисло pi, на которых подход даёт больше "
                f"`floor((n-1)^2/4)` (для L3/L4 — число pi без нужного разреза; "
                f"для L5/L6 — число pi, где достаточное условие не выполнено).\n\n")
        f.write("\n".join(lines) + "\n\nПервые контрпримеры:\n\n" + "\n".join(ex) + "\n\nСемейства (n, минимальный запас):\n\n" + "\n".join(str(r) for r in fam) + "\n")
    print("\n".join(lines))
    print(f"report: {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
