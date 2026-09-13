"""Scan of hypothesis H10 (C27): min_c [2F_c - S_c + R_c] <= B_n + 1.

Version h10_scan-1.0. Construction strict_upper-1.0 (constructions/strict_upper.py),
core oracle-1.0. No distance tables are used; the reference is B_n.

Only the *bound* 2F_c - S_c + R_c is evaluated (cycles, loads, carriers and the
carrier route of strict_upper), not the words themselves, so a full scan of
n = 10 fits in minutes. In N1 the word length equals 2F_c - S_c + H(c) exactly
(C18v) and the carrier-route word has length exactly 2F_c - S_c + R_c
(loss_map, C27v); `--crosscheck` re-verifies that identity here on all
permutations 4 <= n <= 7 by building the words. The second part of H10 (free
reduction, <= B_n) needs the words and is NOT scanned here.

Per permutation pi the scan records, over all c in Z_n:
  min_c   = min_c bound - B_n              (H10, first part)
  mean_c  = mean_c bound - B_n (stored as n*mean - n*B_n, an integer)
  medF    = min over c in argmin_c F_c of bound - B_n   (the "median shift" rule)
  small   = min over c in {0, 1, -1} of bound - B_n
  argmin_in_argminF: whether some shift with minimal bound also minimises F_c
Inputs with min_c >= 0 ("extremal") are stored with their full per-c tables
(rule 15 of AGENTS.md: F_c, S_c, sum M, E_c, theta_c, H(c), R_c, carriers).

Usage:
  python3 experiments/h10_scan.py --crosscheck 7          # words vs bound, 4<=n<=7
  python3 experiments/h10_scan.py --perms 10 --workers 4  # all permutations of n
  python3 experiments/h10_scan.py --families 100          # reflections etc., 4<=n<=N
Output: data/runs/h10_scan/perms_n<N>.json / .md, families.json / .md
"""

import argparse
import collections
import itertools
import json
import os
import random
import sys
import time
from fractions import Fraction
from multiprocessing import Pool

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
for sub in ("oracle", "constructions"):
    sys.path.insert(0, os.path.join(ROOT, sub))

from moves import sigma, rev, CORE_VERSION  # noqa: E402
import strict_upper as su  # noqa: E402

EXPERIMENT_VERSION = "h10_scan-1.0"
OUT_DIR = os.path.join(ROOT, "data", "runs", "h10_scan")


def B(n):
    return n * (n - 1) // 2


def bound_R(pi, c):
    """2F_c - S_c + R_c for the carrier-route variant, with rule-15 statistics."""
    n = len(pi)
    f = [(pi[i] + c) % n for i in range(n)]
    F = S = sumM = E = theta = 0
    carriers = []
    for cyc in su.cycles_of(f):
        cd = su.cycle_data(cyc, n)
        j0 = su.choose_carrier(cyc, cd, n)
        F += cd["F"]
        S += cd["S"]
        sumM += cd["M"]
        E += cd["E"]
        theta += 1 if cd["antipodal_transposition"] else 0
        carriers.append(cyc[j0])
    _, letters = su.carrier_route(carriers, c, n)
    R = len(letters)
    st = {"c": c, "K": len(carriers), "F": F, "S": S, "sumM": sumM, "E": E,
          "theta": theta, "H": su.H_of(c, n), "R": R, "carriers": carriers,
          "bound": 2 * F - S + R}
    return st


def analyze(pi):
    n = len(pi)
    per_c = [bound_R(pi, c) for c in range(n)]
    vals = [st["bound"] for st in per_c]
    Fs = [st["F"] for st in per_c]
    Fmin = min(Fs)
    medc = [c for c in range(n) if Fs[c] == Fmin]
    vmin = min(vals)
    Bn = B(n)
    return {
        "pi": list(pi),
        "min_c": vmin - Bn,
        "argmin_c": vals.index(vmin),
        "n_mean_minus_nB": sum(vals) - n * Bn,
        "medF": min(vals[c] for c in medc) - Bn,
        "small": min(vals[c] for c in (0, 1, n - 1)) - Bn,
        "argmin_in_argminF": any(vals[c] == vmin for c in medc),
        "per_c": per_c,
    }


# ------------------------------------------------------------------ exhaustive scan

def _chunk(args):
    n, first = args
    hist = {k: collections.Counter() for k in ("min_c", "n_mean_minus_nB", "medF", "small")}
    agree = 0
    count = 0
    extremal = []
    worst_mean = []
    rest = [i for i in range(n) if i != first]
    for tail in itertools.permutations(rest):
        pi = (first,) + tail
        a = analyze(pi)
        count += 1
        for k in hist:
            hist[k][a[k]] += 1
        agree += a["argmin_in_argminF"]
        if a["min_c"] >= 0:
            extremal.append(a)
        worst_mean.append((a["n_mean_minus_nB"], a["pi"]))
        if len(worst_mean) > 64:
            worst_mean.sort(reverse=True)
            del worst_mean[8:]
    worst_mean.sort(reverse=True)
    return {"hist": {k: dict(v) for k, v in hist.items()}, "agree": agree,
            "count": count, "extremal": extremal, "worst_mean": worst_mean[:8]}


def scan_perms(n, workers):
    t0 = time.time()
    with Pool(workers) as pool:
        parts = pool.map(_chunk, [(n, first) for first in range(n)])
    hist = {k: collections.Counter() for k in ("min_c", "n_mean_minus_nB", "medF", "small")}
    agree = count = 0
    extremal, worst_mean = [], []
    for p in parts:
        for k in hist:
            for key, v in p["hist"][k].items():
                hist[k][int(key)] += v
        agree += p["agree"]
        count += p["count"]
        extremal += p["extremal"]
        worst_mean += p["worst_mean"]
    worst_mean.sort(reverse=True)
    elapsed = time.time() - t0
    res = {
        "meta": meta(f"--perms {n} --workers {workers}"),
        "n": n, "B_n": B(n), "permutations": count, "shifts": count * n,
        "elapsed_s": round(elapsed, 1), "workers": workers,
        "hist_min_c": sorted(hist["min_c"].items()),
        "hist_medF": sorted(hist["medF"].items()),
        "hist_small": sorted(hist["small"].items()),
        "hist_n_mean_minus_nB": sorted(hist["n_mean_minus_nB"].items()),
        "max_min_c": max(hist["min_c"]),
        "max_mean_c": str(Fraction(max(hist["n_mean_minus_nB"]), n)),
        "argmin_in_argminF": agree,
        "extremal": sorted(extremal, key=lambda a: (-a["min_c"], a["pi"])),
        "worst_mean": worst_mean[:8],
    }
    return res


# ------------------------------------------------------------------ families

def families(n, rng):
    fam = {"sigma_n": tuple(sigma(n)), "rev_n": tuple(rev(n))}
    for h in range(n):
        fam[f"refl_h{h}"] = tuple((h - i) % n for i in range(n))
    for k in range(1, n):
        fam[f"rho_{k}"] = tuple((i + k) % n for i in range(n))
    for k in range(3):
        p = list(range(n))
        rng.shuffle(p)
        fam[f"random_{k}"] = tuple(p)
    for k in range(3):
        p = [(1 - i) % n for i in range(n)]
        for _ in range(k + 1):
            a, b = rng.sample(range(n), 2)
            p[a], p[b] = p[b], p[a]
        fam[f"sigma_perturbed_{k + 1}"] = tuple(p)
    return fam


def scan_families(nmax):
    out = {"meta": meta(f"--families {nmax}"), "by_n": {}}
    for n in range(4, nmax + 1):
        rng = random.Random(20260912 + n)
        rows = []
        for name, pi in families(n, rng).items():
            a = analyze(pi)
            rows.append({"name": name, "pi": list(pi), "min_c": a["min_c"],
                         "argmin_c": a["argmin_c"],
                         "mean_c": str(Fraction(a["n_mean_minus_nB"], n)),
                         "mean_c_float": a["n_mean_minus_nB"] / n,
                         "s": str(Fraction(sum(st["S"] for st in a["per_c"]), n)),
                         "mean_R": str(Fraction(sum(st["R"] for st in a["per_c"]), n)),
                         "mean_F": str(Fraction(sum(st["F"] for st in a["per_c"]), n)),
                         "medF": a["medF"], "small": a["small"]})
        worst_mean = max(rows, key=lambda r: r["mean_c_float"])
        out["by_n"][n] = {
            "B_n": B(n), "rows": rows,
            "max_min_c": max(r["min_c"] for r in rows),
            "argmax_min_c": [r["name"] for r in rows if r["min_c"] == max(x["min_c"] for x in rows)],
            "max_mean_c": worst_mean["mean_c"], "argmax_mean_c": worst_mean["name"],
        }
        print(f"n={n}: max min_c={out['by_n'][n]['max_min_c']} ({out['by_n'][n]['argmax_min_c'][:3]}), "
              f"max mean_c={worst_mean['mean_c']} ({worst_mean['name']})", flush=True)
    return out


# ------------------------------------------------------------------ crosscheck against words

def crosscheck(nmax):
    """bound_R == len(word) of strict_upper.shift_word(route='carrier'), all pi, all c."""
    total = 0
    for n in range(4, nmax + 1):
        for pi in itertools.permutations(range(n)):
            for c in range(n):
                sw = su.shift_word(pi, c, "carrier")
                st = bound_R(pi, c)
                assert sw["stats"]["len"] == st["bound"], (pi, c)
                assert sw["stats"]["F_c"] == st["F"] and sw["stats"]["S_c"] == st["S"]
                assert sw["stats"]["route_len"] == st["R"]
                assert sorted(sw["stats"]["carriers"]) == sorted(st["carriers"])
                total += 1
        print(f"crosscheck n={n} OK", flush=True)
    return total


# ------------------------------------------------------------------ output

def meta(args):
    return {"experiment": EXPERIMENT_VERSION, "construction": su.CONSTRUCTION_VERSION,
            "core": CORE_VERSION, "command": f"python3 experiments/h10_scan.py {args}",
            "python": sys.version.split()[0]}


def table_md(a, n):
    lines = ["| c | K | F_c | S_c | sum M | E_c | theta | H | R | carriers | 2F-S+R | -B_n |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for st in a["per_c"]:
        lines.append(f"| {st['c']} | {st['K']} | {st['F']} | {st['S']} | {st['sumM']} | {st['E']} | "
                     f"{st['theta']} | {st['H']} | {st['R']} | {','.join(map(str, st['carriers']))} | "
                     f"{st['bound']} | {st['bound'] - B(n)} |")
    return "\n".join(lines)


def write_perms(res):
    n = res["n"]
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(os.path.join(OUT_DIR, f"perms_n{n}.json"), "w") as f:
        json.dump(res, f, indent=1)
    md = [f"# H10 scan, all permutations, n = {n}", "",
          f"{res['meta']['command']}; {res['meta']['construction']}, {res['meta']['core']}; "
          f"{res['permutations']} permutations, {res['shifts']} shifts, {res['elapsed_s']} s, "
          f"{res['workers']} workers. B_n = {res['B_n']}.", "",
          f"- `min_c [2F_c - S_c + R_c] - B_n`: {res['hist_min_c']}; max = {res['max_min_c']}",
          f"- median-F rule (`c in argmin F_c`): {res['hist_medF']}",
          f"- `c in {{0, 1, -1}}`: {res['hist_small']}",
          f"- `mean_c - B_n` max = {res['max_mean_c']}; some minimiser also minimises F_c on "
          f"{res['argmin_in_argminF']} of {res['permutations']}",
          "", f"## Extremal inputs (min_c >= 0): {len(res['extremal'])}", ""]
    for a in res["extremal"]:
        md += [f"### pi = {tuple(a['pi'])}, min_c - B_n = {a['min_c']} at c = {a['argmin_c']}", "",
               table_md(a, n), ""]
    md += ["## Largest mean_c (n*mean - n*B_n, pi)", ""]
    md += [f"- {v}: {tuple(p)}" for v, p in res["worst_mean"]]
    with open(os.path.join(OUT_DIR, f"perms_n{n}.md"), "w") as f:
        f.write("\n".join(md) + "\n")


def write_families(out):
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(os.path.join(OUT_DIR, "families.json"), "w") as f:
        json.dump(out, f, indent=1)
    md = ["# H10 scan, families", "", out["meta"]["command"], "",
          "| n | B_n | max min_c | where | max mean_c | where | sigma_n: min_c, mean_c, s, mean R |",
          "|---|---|---|---|---|---|---|"]
    for n, r in out["by_n"].items():
        sg = next(x for x in r["rows"] if x["name"] == "sigma_n")
        md.append(f"| {n} | {r['B_n']} | {r['max_min_c']} | {','.join(r['argmax_min_c'][:4])} | "
                  f"{r['max_mean_c']} | {r['argmax_mean_c']} | {sg['min_c']}, {sg['mean_c']}, {sg['s']}, {sg['mean_R']} |")
    with open(os.path.join(OUT_DIR, "families.md"), "w") as f:
        f.write("\n".join(md) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perms", type=int)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--families", type=int)
    ap.add_argument("--crosscheck", type=int)
    a = ap.parse_args()
    if a.crosscheck:
        t = crosscheck(a.crosscheck)
        print(f"crosscheck: {t} (pi, c) pairs, 4<=n<={a.crosscheck}: word length == 2F_c - S_c + R_c")
    if a.families:
        write_families(scan_families(a.families))
    if a.perms:
        res = scan_perms(a.perms, a.workers)
        write_perms(res)
        print(f"n={a.perms}: max min_c - B_n = {res['max_min_c']}, hist {res['hist_min_c']}, "
              f"max mean_c - B_n = {res['max_mean_c']}, {res['elapsed_s']} s")


if __name__ == "__main__":
    main()
