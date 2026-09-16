#!/usr/bin/env python3
"""h13i_probes.py — the auxiliary tables of docs/notes/h13i_attempt.md.

Reproduces, for the exhaustive ranges stated in the note:
  arc   — how many toric class representatives satisfy the arc criterion (S),
          and max I among those that do not (section 4);
  fam   — the n^2 rigid families F(m, t) of n cuts each (section 3.3);
  lines — averaging inv over a single line of the cut torus (section 3.1);
  dg    — the Diaconis-Graham route min (D - T) and the plain footrule min D
          (section 3.4);
  affine— I for pi(i) = c*i + d as a function of c (section 5).

Usage: python3 experiments/h13i_probes.py [--nmax 8] [--out probes.json]
Version h13i_probes-1.0.  Independent of the C scanner; plain O(n^2) inversion
counting everywhere except the O(1)-update table, which check_H13I.py verifies.
"""
import argparse
import itertools
import json
from math import gcd


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def line(pi, n, a, b):
    w = [0] * n
    for i in range(n):
        w[(i - a) % n] = (pi[i] - b) % n
    return w


def inv_table(pi, n):
    """inv(a, b) for all cuts, via the shift identities (lemma 1.5)."""
    pinv = [0] * n
    for i, v in enumerate(pi):
        pinv[v] = i
    tab = [[0] * n for _ in range(n)]
    cur_a0 = inversions(list(pi))
    for a in range(n):
        cur = cur_a0
        for b in range(n):
            tab[a][b] = cur
            cur += (n - 1) - 2 * ((pinv[b] - a) % n)
        cur_a0 += (n - 1) - 2 * pi[a]
    return tab


def bound(n):
    return (n - 1) ** 2 // 4


def reps(n):
    for rest in itertools.permutations(range(1, n)):
        yield (0,) + rest


def has_arc(pi, n, k):
    for a in range(n):
        s = set(pi[(a + t) % n] for t in range(k))
        for b in range(n):
            if all(((b + t) % n) in s for t in range(k)):
                return True
    return False


def cayley(w):
    n = len(w)
    seen, c = [False] * n, 0
    for i in range(n):
        if not seen[i]:
            c += 1
            j = i
            while not seen[j]:
                seen[j] = True
                j = w[j]
    return n - c


def probe_arc(n):
    k, with_s, max_no, arg_no = n // 2, 0, -1, None
    for pi in reps(n):
        I = min(min(r) for r in inv_table(pi, n))
        if has_arc(pi, n, k):
            with_s += 1
        elif I > max_no:
            max_no, arg_no = I, pi
    return {"k": k, "with_S": with_s, "max_I_without_S": max_no,
            "arg_without_S": arg_no}


def probe_fam(n):
    worst = {(m, t): -1 for m in range(n) for t in range(n)}
    for pi in reps(n):
        tab = inv_table(pi, n)
        for (m, t) in worst:
            v = min(tab[a][(pi[(a + m) % n] - t) % n] for a in range(n))
            if v > worst[(m, t)]:
                worst[(m, t)] = v
    best = min(worst.values())
    return {"min_over_families_of_max": best, "sufficient_families":
            [list(f) for f in worst if worst[f] <= bound(n)]}


def probe_lines(n):
    worst, arg = -1, None
    for pi in reps(n):
        tab = inv_table(pi, n)
        best = min(sum(tab[a][(a * g + t) % n] for a in range(n))
                   for g in range(n) for t in range(n))
        best = min(best, min(sum(tab[t][b] for b in range(n)) for t in range(n)))
        if best > worst:
            worst, arg = best, pi
    return {"max_over_pi_of_min_line_sum": worst, "as_average": worst / n,
            "arg": arg}


def probe_dg(n):
    wd, wdt, argdt = -1, -1, None
    for pi in reps(n):
        bd, bdt = 10 ** 9, 10 ** 9
        for a in range(n):
            for b in range(n):
                w = line(pi, n, a, b)
                D = sum(abs(j - w[j]) for j in range(n))
                bd = min(bd, D)
                bdt = min(bdt, D - cayley(w))
        wd = max(wd, bd)
        if bdt > wdt:
            wdt, argdt = bdt, pi
    return {"max_min_D": wd, "max_min_D_minus_T": wdt, "arg": argdt}


def probe_affine(n):
    out = {}
    for c in range(1, n):
        if gcd(c, n) != 1:
            continue
        out[c] = min(inversions([(c * x + e) % n for x in range(n)])
                     for e in range(n))
    others = [v for c, v in out.items() if c not in (1, n - 1)]
    return {"I_by_c": out, "max_I": max(out.values()),
            "max_I_c_not_pm1": max(others) if others else None,
            "slack_c_not_pm1": (bound(n) - max(others)) if others else None}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--affine-ns", default="16,17,18,20,24,30,31,32,40,50,60,64,101,128")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    res = {}
    for n in range(4, args.nmax + 1):
        row = {"bound": bound(n), "arc": probe_arc(n), "fam": probe_fam(n),
               "lines": probe_lines(n), "dg": probe_dg(n)}
        res[str(n)] = row
        print(f"n={n} bound={row['bound']} arc={row['arc']} "
              f"fam_min={row['fam']['min_over_families_of_max']} "
              f"line_avg={row['lines']['as_average']:.3f} "
              f"minD={row['dg']['max_min_D']} minD-T={row['dg']['max_min_D_minus_T']}")
    aff = {}
    for n in (int(x) for x in args.affine_ns.split(",") if x.strip()):
        aff[str(n)] = probe_affine(n)
        print(f"affine n={n} bound={bound(n)} max_I={aff[str(n)]['max_I']} "
              f"max_I(c!=+-1)={aff[str(n)]['max_I_c_not_pm1']} "
              f"slack={aff[str(n)]['slack_c_not_pm1']}")
    res["affine"] = aff
    if args.out:
        with open(args.out, "w") as f:
            json.dump(res, f, indent=1)


if __name__ == "__main__":
    main()
