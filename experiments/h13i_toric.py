"""H13-I (toric reformulation of the double cut): pairwise inversion-count
identity and three weakened cut families, tested against I(pi) and the
target floor((n-1)^2/4).

Definitions (docs/notes/h13_line_model.md §0): for pi in S_n (as a bijection
Z_n -> Z_n) and a cut (a, b) in Z_n x Z_n, the line is
w_j = (pi((j + a) mod n) - b) mod n, j = 0..n-1; I(pi) = min_{a,b} inv(w).

This script:

  1. Verifies, by direct computation, the identity

       sum_{a,b in Z_n} inv(w^{a,b}) = sum_{i<i'} [ n(d + e) - 2 d e ]

     where d = i' - i, e = (pi(i') - pi(i)) mod n (unordered pairs i<i' as
     integers 0 <= i < i' <= n-1).  Proof: for a fixed pair {i, i'}, the
     number of a in Z_n with i preceding i' in the rotated reading is n - d,
     and with i' preceding i is d (standard cyclic-cut count); independently,
     the number of b with (pi(i)-b) mod n > (pi(i')-b) mod n is e, and the
     complementary count is n - e.  The pair is inverted for (a, b) exactly
     when (i before i' and value(i) > value(i')) or (i' before i and
     value(i') > value(i)), giving (n-d) e + d (n-e) = n(d+e) - 2de choices
     of (a, b) out of n^2.  Summing over pairs gives the identity.  This is
     an exact algebraic identity, verified here only as a numerical sanity
     check (assert), not as a proof technique for H13-I: the averaging
     min <= (sum over a,b)/n^2 that it would support is exactly the
     averaging approach already noted as insufficient (docs/notes/
     h13_line_model.md, PLAN §8) and is not used as a bound below.

  2. Tests three weaker cut families as *candidate sufficient schemes* for
     H13-I (each would have given a much shorter proof if it worked):

       single_a:   b = 0 fixed, only a varies (n candidates);
       a_then_avg: a chosen to minimize J(h_a) = sum_{j<k} (h_a(k)-h_a(j))
                   mod n (the sum over b of inv(h_a shifted by b), so
                   floor(J(h_a)/n) upper-bounds min_b inv(h_a, b) by
                   pigeonhole), then that bound is compared to the target;
       vertex_cut: (a, b) = (i0, pi(i0)) for some data point i0 (n
                   candidates; places that point at value 0, contributing no
                   inversions), best of the n choices per pi.

     All three are REFUTED as sufficient (exhaustive smallest failing pi
     reported); genuine two-parameter, non-vertex cuts are required, as the
     project's own H13-I attempts already found for averaging, one-element
     induction, and the "first element = t" family
     (docs/notes/h13_line_model.md §1.7, §6).

Usage: python3 experiments/h13i_toric.py [--amax 8] [--formula-trials 200]
Output: data/runs/h13i_toric/report.json, report.md
"""

import argparse
import itertools
import json
import math
import os
import random
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "data", "runs", "h13i_toric")

VERSION = "h13i_toric-1.0"


def inv_count(seq):
    n = len(seq)
    c = 0
    for j in range(n):
        for k in range(j + 1, n):
            if seq[j] > seq[k]:
                c += 1
    return c


def line(f, a, b):
    n = len(f)
    return [(f[(j + a) % n] - b) % n for j in range(n)]


def I_of_pi(f):
    n = len(f)
    best = None
    for a in range(n):
        for b in range(n):
            iv = inv_count(line(f, a, b))
            if best is None or iv < best:
                best = iv
    return best


def total_inv_over_ab(f):
    n = len(f)
    return sum(inv_count(line(f, a, b)) for a in range(n) for b in range(n))


def formula_sum(f):
    n = len(f)
    s = 0
    for i in range(n):
        for ip in range(i + 1, n):
            d = ip - i
            e = (f[ip] - f[i]) % n
            s += n * (d + e) - 2 * d * e
    return s


def single_a_best(f):
    n = len(f)
    return min(inv_count(line(f, a, 0)) for a in range(n))


def J_of_ha(f, a):
    n = len(f)
    h = [f[(j + a) % n] for j in range(n)]
    s = 0
    for j in range(n):
        for k in range(j + 1, n):
            s += (h[k] - h[j]) % n
    return s


def a_then_avg_bound(f):
    n = len(f)
    minJ = min(J_of_ha(f, a) for a in range(n))
    return minJ // n


def vertex_cut_best(f):
    n = len(f)
    return min(inv_count(line(f, i0, f[i0])) for i0 in range(n))


def check_formula(nmax, trials, rng):
    results = []
    for n in range(3, nmax + 1):
        ok = True
        checked = 0
        for _ in range(trials):
            f = list(range(n))
            rng.shuffle(f)
            t1 = total_inv_over_ab(f)
            t2 = formula_sum(f)
            checked += 1
            if t1 != t2:
                ok = False
                results.append({"n": n, "ok": False, "pi": f, "direct": t1, "formula": t2})
                break
        if ok:
            results.append({"n": n, "ok": True, "checked": checked})
    return results


def check_family(name, fn, nmax):
    results = []
    for n in range(3, nmax + 1):
        target = ((n - 1) ** 2) // 4
        worst = -1
        worst_pi = None
        for perm in itertools.permutations(range(n)):
            f = list(perm)
            v = fn(f)
            if v > worst:
                worst = v
                worst_pi = f
        results.append({
            "n": n,
            "target": target,
            "max_over_pi": worst,
            "sufficient": worst <= target,
            "worst_pi": worst_pi if worst > target else None,
        })
        if worst > target:
            # smallest failing n found; stop (rule 6: minimal counterexample)
            break
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amax", type=int, default=8, help="max n for exhaustive family checks")
    ap.add_argument("--formula-nmax", type=int, default=9)
    ap.add_argument("--formula-trials", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()

    formula = check_formula(args.formula_nmax, args.formula_trials, rng)
    assert all(r["ok"] for r in formula), "pairwise inversion-count identity failed"

    families = {
        "single_a (b=0 fixed)": check_family("single_a", single_a_best, args.amax),
        "a_then_avg (best a, pigeonhole on b)": check_family("a_then_avg", a_then_avg_bound, args.amax),
        "vertex_cut (best data point as origin)": check_family("vertex_cut", vertex_cut_best, args.amax),
    }

    elapsed = time.time() - t0
    report = {
        "version": VERSION,
        "elapsed_s": elapsed,
        "formula_identity": formula,
        "families": families,
    }
    with open(os.path.join(OUT, "report.json"), "w") as fh:
        json.dump(report, fh, indent=2)

    with open(os.path.join(OUT, "report.md"), "w") as fh:
        fh.write("# H13-I toric reformulation: identity check + three refuted families\n\n")
        fh.write(f"Version {VERSION}, elapsed {elapsed:.1f}s.\n\n")
        fh.write("## Pairwise inversion-count identity\n\n")
        fh.write("`sum_{a,b} inv(w^{a,b}) = sum_{i<i'} [n(d+e) - 2de]` verified by direct\n")
        fh.write("computation (random permutations, `--formula-trials` each) for all n in\n")
        fh.write(f"range checked (3..{args.formula_nmax}): all OK.\n\n")
        fh.write("## Weakened cut families (exhaustive over all pi, per n, up to the smallest\n")
        fh.write("failing n or `--amax`)\n\n")
        for name, rows in families.items():
            fh.write(f"### {name}\n\n")
            fh.write("| n | target floor((n-1)^2/4) | max over pi | sufficient? | worst pi |\n")
            fh.write("|---|---|---|---|---|\n")
            for r in rows:
                fh.write(f"| {r['n']} | {r['target']} | {r['max_over_pi']} | {r['sufficient']} | {r['worst_pi']} |\n")
            fh.write("\n")
    print(json.dumps({"elapsed_s": elapsed, "formula_ok": all(r["ok"] for r in formula)}, indent=2))


if __name__ == "__main__":
    main()
