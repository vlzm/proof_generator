"""Toric correlation identity for the cut-inversion function (H13-I).

For pi on Z_n, a cut (q, c) gives the line w_j = (pi(q+1+j) - (q+1-c)) mod n
(PLAN, H13-I).  Put s = (q - c) mod n and

    sig(t) = n - 1 - 2 * (t mod n)          (centred sawtooth, sum 0)
    S(q,s) = sum_i sig(q - i) * sig(s - pi(i))
    E(pi)  = sum_{i,j} sig(j - i) * sig(pi(j) - pi(i))

Theorem A (identity):  4 n^2 inv(q,c) = n^2 (n^2 - 1) - E - 2 n S(q,s).
Theorem B (energy):    E >= 2 n^2 (n-1) - n^2 (n^2-1)/3, equality iff pi is a
                       reflection pi_h(i) = h - i.
Corollary C:           avg_{q,c} inv = (n^2-1)/4 - E/(4n^2) <= (n-1)(2n-1)/6,
                       hence I(pi) <= floor((n-1)(2n-1)/6) for every pi.

This script checks A, B, C exhaustively and then tests the candidate sufficient
conditions for H13-I that the identity suggests (all of them fail; see
docs/notes/h13i_verdict.md).  Version toric_corr-1.0.

Usage: python3 experiments/toric_correlation.py [--nmax 8] [--idmax 7]
Output: data/runs/toric_correlation/report.json, report.md
"""

import argparse
import itertools
import json
import os
import random
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUTDIR = os.path.join(ROOT, "data", "runs", "toric_correlation")
VERSION = "toric_corr-1.0"


# ---------------------------------------------------------------- primitives

def line(perm, q, c):
    """The H13-I line of the cut (q, c)."""
    n = len(perm)
    return [(perm[(q + 1 + j) % n] - (q + 1 - c)) % n for j in range(n)]


def inv_count(w):
    n = len(w)
    return sum(1 for j in range(n) for k in range(j + 1, n) if w[j] > w[k])


def s_table(perm):
    """S(q, s) for all q, s in Z_n (O(n^3), plain)."""
    n = len(perm)
    sig = [n - 1 - 2 * t for t in range(n)]
    out = [[0] * n for _ in range(n)]
    for i, v in enumerate(perm):
        a = [sig[(q - i) % n] for q in range(n)]
        b = [sig[(s - v) % n] for s in range(n)]
        for q in range(n):
            aq = a[q]
            row = out[q]
            for s in range(n):
                row[s] += aq * b[s]
    return out


def energy(perm):
    n = len(perm)
    sig = [n - 1 - 2 * t for t in range(n)]
    return sum(sig[(j - i) % n] * sig[(perm[j] - perm[i]) % n]
               for i in range(n) for j in range(n))


def energy_bound(n):
    """2 n^2 (n-1) - n^2 (n^2-1)/3 (an integer for every n)."""
    assert (n * n * (n * n - 1)) % 3 == 0
    return 2 * n * n * (n - 1) - n * n * (n * n - 1) // 3


def avg_bound(n):
    """floor((n-1)(2n-1)/6) -- the proved upper bound on I(pi)."""
    return ((n - 1) * (2 * n - 1)) // 6


def reflections(n):
    return set(tuple((h - i) % n for i in range(n)) for h in range(n))


def peak_of_reflection(n):
    """max_{q,s} S for a reflection: -n(n^2-1)/3 + 2 n floor(n^2/4)."""
    return -n * (n * n - 1) // 3 + 2 * n * (n * n // 4)


# ---------------------------------------------------------------- the parts

def part_identity(idmax):
    """Theorem A against brute force, all pi and all n^2 cuts."""
    res = {}
    for n in range(2, idmax + 1):
        bad = 0
        checked = 0
        for perm in itertools.permutations(range(n)):
            E = energy(perm)
            S = s_table(perm)
            for q in range(n):
                for c in range(n):
                    inv = inv_count(line(perm, q, c))
                    s = (q - c) % n
                    if 4 * n * n * inv != n * n * (n * n - 1) - E - 2 * n * S[q][s]:
                        bad += 1
                    checked += 1
        res[n] = {"cuts_checked": checked, "failures": bad}
    return res


def part_energy(nmax):
    """Theorem B: the bound and its equality set."""
    res = {}
    for n in range(2, nmax + 1):
        bnd = energy_bound(n)
        below = 0
        eq = set()
        emin = None
        for perm in itertools.permutations(range(n)):
            E = energy(perm)
            if emin is None or E < emin:
                emin = E
            if E < bnd:
                below += 1
            elif E == bnd:
                eq.add(perm)
        res[n] = {"bound": bnd, "min_E": emin, "below_bound": below,
                  "equality_count": len(eq),
                  "equality_is_reflections": eq == reflections(n)}
    return res


def part_bound(nmax):
    """Corollary C against the exhaustive I(pi) and max_pi I(pi)."""
    res = {}
    for n in range(4, nmax + 1):
        bnd = avg_bound(n)
        target = (n - 1) ** 2 // 4
        max_I = 0
        max_avg_num = None          # 4 n^2 * avg inv, integer
        viol = 0
        for perm in itertools.permutations(range(n)):
            E = energy(perm)
            S = s_table(perm)
            smax = max(max(row) for row in S)
            # inv_min = ((n^2-1) n^2 - E - 2 n smax) / (4 n^2)
            num = n * n * (n * n - 1) - E - 2 * n * smax
            assert num % (4 * n * n) == 0
            I = num // (4 * n * n)
            max_I = max(max_I, I)
            avg_num = n * n * (n * n - 1) - E
            if max_avg_num is None or avg_num > max_avg_num:
                max_avg_num = avg_num
            if I > bnd:
                viol += 1
        res[n] = {"proved_bound_floor": bnd, "max_I": max_I,
                  "target_floor_(n-1)^2/4": target,
                  "max_avg_inv_times_4n2": max_avg_num,
                  "avg_bound_times_4n2": 4 * n * n * ((n - 1) * (2 * n - 1)) // 6,
                  "violations": viol}
    return res


def part_candidates(nmax):
    """The candidate sufficient conditions suggested by the identity.

    (P) peak bound            max S >= peak_of_reflection(n)
    (M) 2 max S + min S >= 2 thr
    (T) three-point (n even)  sum_j [S(j,pi(j)) + S(j+n/2,pi(j)) + S(j,pi(j)+n/2)] >= 2 n thr
    (L) slope family          max over lam, s of (1/n) sum_a S(a, lam a + s)
    (W) shifted-graph family  max over (s,t) of (1/n) sum_j S(j+s, pi(j)+t)
    where thr = n floor(n/2) + n(n-1)/2 and the exact criterion is
    2 n max S + E >= 2 n thr.
    """
    res = {}
    for n in range(4, nmax + 1):
        thr = n * (n // 2) + n * (n - 1) // 2
        peak = peak_of_reflection(n)
        cnt = {"P": 0, "M": 0, "T": 0, "L": 0, "W": 0, "exact": 0}
        first = {}
        total = 0
        # precompute slope kernels D_lam(m) = sum_a sig(a) sig(lam a + m)
        sig = [n - 1 - 2 * t for t in range(n)]
        D = [[sum(sig[a] * sig[(lam * a + m) % n] for a in range(n))
              for m in range(n)] for lam in range(n)]
        for perm in itertools.permutations(range(n)):
            total += 1
            E = energy(perm)
            S = s_table(perm)
            smax = max(max(r) for r in S)
            smin = min(min(r) for r in S)
            if 2 * n * smax + E < 2 * n * thr:
                cnt["exact"] += 1
                first.setdefault("exact", list(perm))
            if smax < peak:
                cnt["P"] += 1
                first.setdefault("P", list(perm))
            if 2 * smax + smin < 2 * thr:
                cnt["M"] += 1
                first.setdefault("M", list(perm))
            if n % 2 == 0:
                h = n // 2
                A = sum(S[j][perm[j]] + S[(j + h) % n][perm[j]]
                        + S[j][(perm[j] + h) % n] for j in range(n))
                if A < 2 * n * thr:
                    cnt["T"] += 1
                    first.setdefault("T", list(perm))
            best_l = max(max(sum(D[lam][(lam * i - perm[i] + s) % n]
                                 for i in range(n)) for s in range(n))
                         for lam in range(n))
            if 2 * best_l + E < 2 * n * thr:
                cnt["L"] += 1
                first.setdefault("L", list(perm))
            best_w = max(max(sum(S[(j + sh) % n][(perm[j] + t) % n]
                                 for j in range(n)) for t in range(n))
                         for sh in range(n))
            if 2 * best_w + E < 2 * n * thr:
                cnt["W"] += 1
                first.setdefault("W", list(perm))
        res[n] = {"perms": total, "thr": thr, "reflection_peak": peak,
                  "failures": cnt, "first_failure": first}
    return res


def part_random_peak(ns, trials, seed):
    """The peak bound on random permutations at larger n."""
    rng = random.Random(seed)
    out = {}
    for n in ns:
        need = peak_of_reflection(n)
        worst = None
        worst_perm = None
        for _ in range(trials):
            p = list(range(n))
            rng.shuffle(p)
            smax = max(max(r) for r in s_table(p))
            if worst is None or smax < worst:
                worst, worst_perm = smax, list(p)
        out[n] = {"reflection_peak": need, "min_random_maxS": worst,
                  "fails": worst < need, "witness": worst_perm if worst < need else None}
    return out


# ---------------------------------------------------------------- driver

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--idmax", type=int, default=6)
    ap.add_argument("--random-n", type=int, nargs="*", default=[12, 20, 30, 50])
    ap.add_argument("--trials", type=int, default=30)
    ap.add_argument("--seed", type=int, default=20260913)
    args = ap.parse_args()

    t0 = time.time()
    report = {"version": VERSION, "args": vars(args)}
    report["identity"] = part_identity(args.idmax)
    report["energy"] = part_energy(min(args.nmax, 8))
    report["bound"] = part_bound(min(args.nmax, 8))
    report["candidates"] = part_candidates(min(args.nmax, 8))
    report["random_peak"] = part_random_peak(args.random_n, args.trials, args.seed)
    report["seconds"] = round(time.time() - t0, 1)

    os.makedirs(OUTDIR, exist_ok=True)
    with open(os.path.join(OUTDIR, "report.json"), "w") as f:
        json.dump(report, f, indent=1, sort_keys=True)

    lines = ["# Toric correlation identity — report", "",
             f"Version {VERSION}; {report['seconds']} s.", "",
             "## Theorem A (identity), brute force over all cuts", ""]
    for n, d in sorted(report["identity"].items()):
        lines.append(f"- n = {n}: {d['cuts_checked']} cuts, failures {d['failures']}")
    lines += ["", "## Theorem B (energy bound) and its equality set", ""]
    for n, d in sorted(report["energy"].items()):
        lines.append(f"- n = {n}: bound {d['bound']}, min E {d['min_E']}, "
                     f"below {d['below_bound']}, equality on {d['equality_count']} "
                     f"perms, = reflections: {d['equality_is_reflections']}")
    lines += ["", "## Corollary C: I(pi) <= floor((n-1)(2n-1)/6)", ""]
    for n, d in sorted(report["bound"].items()):
        lines.append(f"- n = {n}: proved bound {d['proved_bound_floor']}, "
                     f"max_pi I = {d['max_I']}, target floor((n-1)^2/4) = "
                     f"{d['target_floor_(n-1)^2/4']}, violations {d['violations']}")
    lines += ["", "## Candidate sufficient conditions for H13-I (all fail)", ""]
    for n, d in sorted(report["candidates"].items()):
        lines.append(f"- n = {n} ({d['perms']} perms), thr = {d['thr']}, "
                     f"reflection peak = {d['reflection_peak']}: failures "
                     + ", ".join(f"{k}={v}" for k, v in sorted(d["failures"].items())))
        for k, p in sorted(d["first_failure"].items()):
            lines.append(f"    first {k}-failure: {p}")
    lines += ["", "## Peak bound on random permutations", ""]
    for n, d in sorted(report["random_peak"].items()):
        lines.append(f"- n = {n}: reflection peak {d['reflection_peak']}, "
                     f"min random max S {d['min_random_maxS']}, fails: {d['fails']}")
    with open(os.path.join(OUTDIR, "report.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
