"""toric_sampling-1.0 — structured and random samples for H13-I at larger n.

H13-I:  I(pi) = min over the n^2 double cuts (a,s) of inv(w) <= floor((n-1)^2/4),
        w_j = (pi(a+j) - s) mod n.

Reported for every n: the value of I on the structured families (reflections,
rotations, affine pi(i) = a*i + b, reversals of blocks), the maximum of I over
random samples and over hill-climbing runs that try to maximise I, and whether
any near-extremal permutation is something other than a reflection.

Also reported per permutation: c(pi) = (1/n^2) sum_{i<k} (T(k-i) - T(pi(k)-pi(i)))^2
(toric invariant of the identity inv = Q/n + c, docs/proofs/C37_toric_identity.md),
min_r J(r) (the aligned relaxation) and the row statistic sum_a min_s inv(a,s)
of conjecture P'.

Usage: python3 experiments/toric_sampling.py --ns 12,15,20,30 --samples 200 --seed 1
"""
from __future__ import annotations

import argparse
import json
import math
import random
import sys
from fractions import Fraction


# ---------------------------------------------------------------- core


def inv_linear(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def inv_table(pi):
    """tab[a][s] = inv of the line of the double cut (a, s).

    Filled with the increment recursions
        inv(a+1,s) - inv(a,s) = (n-1) - 2*((pi(a)-s) mod n)
        inv(a,s+1) - inv(a,s) = (n-1) - 2*((pi^{-1}(s)-a) mod n)
    (checked against the definition by checks/check_C37.py).
    """
    n = len(pi)
    ipi = [0] * n
    for i, v in enumerate(pi):
        ipi[v] = i
    tab = [[0] * n for _ in range(n)]
    tab[0][0] = inv_linear(list(pi))
    for s in range(n - 1):
        tab[0][s + 1] = tab[0][s] + (n - 1) - 2 * ipi[s]
    for s in range(n):
        row0 = tab[0][s]
        prev = row0
        for a in range(n - 1):
            prev = prev + (n - 1) - 2 * ((pi[a] - s) % n)
            tab[a + 1][s] = prev
    return tab


def toric_I(pi):
    tab = inv_table(pi)
    return min(min(row) for row in tab)


def c_invariant(pi):
    """c(pi) = (1/n^2) * sum_{i<k} (T(k-i) - T(pi(k)-pi(i)))^2."""
    n = len(pi)
    tot = 0
    for i in range(n):
        for k in range(i + 1, n):
            tot += ((k - i) % n - (pi[k] - pi[i]) % n) ** 2
    return Fraction(tot, n * n)


def min_J(pi):
    """min over r of J(r) = c + (1/n) sum_i D(d_i - r)^2, d_i = (i-pi(i)) mod n."""
    n = len(pi)
    c = c_invariant(pi)
    best = None
    for r in range(n):
        sq = 0
        for i in range(n):
            d = (i - pi[i] - r) % n
            sq += min(d, n - d) ** 2
        val = c + Fraction(sq, n)
        if best is None or val < best:
            best = val
    return best


def row_statistic(pi):
    """sum_a min_s inv(a,s) — conjecture P' says this is <= n*floor((n-1)^2/4)."""
    tab = inv_table(pi)
    return sum(min(row) for row in tab)


# ---------------------------------------------------------------- families


def reflections(n):
    return [tuple((h - i) % n for i in range(n)) for h in range(n)]


def rotations(n):
    return [tuple((i + b) % n for i in range(n)) for b in range(n)]


def affine(n):
    out = []
    for a in range(1, n):
        if math.gcd(a, n) != 1:
            continue
        for b in (0, 1, n // 2, n - 1):
            out.append((a, b, tuple((a * i + b) % n for i in range(n))))
    return out


def reflection_perturbations(n, k, rng, count):
    """reflections with k random transpositions of values applied."""
    out = []
    for _ in range(count):
        h = rng.randrange(n)
        pi = [(h - i) % n for i in range(n)]
        for _ in range(k):
            i, j = rng.sample(range(n), 2)
            pi[i], pi[j] = pi[j], pi[i]
        out.append(tuple(pi))
    return out


def hill_climb(n, rng, restarts, sweeps):
    """greedy search for permutations with large I(pi) (transposition moves)."""
    best = (-1, None)
    for _ in range(restarts):
        pi = list(range(n))
        rng.shuffle(pi)
        cur = toric_I(pi)
        for _ in range(sweeps):
            improved = False
            for i in range(n):
                for j in range(i + 1, n):
                    pi[i], pi[j] = pi[j], pi[i]
                    val = toric_I(pi)
                    if val > cur:
                        cur = val
                        improved = True
                    else:
                        pi[i], pi[j] = pi[j], pi[i]
            if not improved:
                break
        if cur > best[0]:
            best = (cur, tuple(pi))
    return best


def is_reflection(pi):
    n = len(pi)
    h = (pi[0]) % n
    return all((pi[i] + i) % n == h for i in range(n))


# ---------------------------------------------------------------- driver


def run(n, samples, seed, restarts, sweeps):
    rng = random.Random(seed * 1000 + n)
    bound = (n - 1) ** 2 // 4
    rep = {"n": n, "bound": bound, "n_times_bound": n * bound}

    refl = [toric_I(p) for p in reflections(n)]
    rep["reflections_I"] = sorted(set(refl))
    rep["reflections_all_equal_bound"] = all(v == bound for v in refl)
    rep["reflections_c"] = str(c_invariant(reflections(n)[0]))
    rep["c_max_theory"] = str(Fraction((n - 1) * (n - 2), 6))

    rep["rotations_I"] = sorted(set(toric_I(p) for p in rotations(n)))

    aff = [(a, b, toric_I(p), is_reflection(p)) for a, b, p in affine(n)]
    non_refl = [x for x in aff if not x[3]]
    rep["affine_count"] = len(aff)
    rep["affine_max_I_nonreflection"] = max((x[2] for x in non_refl), default=None)
    rep["affine_argmax_nonreflection"] = max(non_refl, key=lambda x: x[2])[:3] if non_refl else None

    worst_rand = (-1, None)
    for _ in range(samples):
        pi = list(range(n))
        rng.shuffle(pi)
        v = toric_I(pi)
        if v > worst_rand[0]:
            worst_rand = (v, tuple(pi))
    rep["random_samples"] = samples
    rep["random_max_I"] = worst_rand[0]

    for k in (1, 2, 3):
        pert = reflection_perturbations(n, k, rng, max(20, samples // 4))
        pert = [p for p in pert if not is_reflection(p)]
        vals = [toric_I(p) for p in pert]
        rep[f"reflection_plus_{k}_transpositions_max_I"] = max(vals) if vals else None
        rep[f"reflection_plus_{k}_transpositions_count"] = len(vals)

    hc = hill_climb(n, rng, restarts, sweeps)
    rep["hillclimb_restarts"] = restarts
    rep["hillclimb_max_I"] = hc[0]
    rep["hillclimb_argmax_is_reflection"] = is_reflection(hc[1]) if hc[1] else None
    rep["hillclimb_argmax"] = list(hc[1]) if hc[1] else None

    # violations of H13-I anywhere in this sample, and of conjecture P'
    pool = [p for p in reflections(n)] + [p for p in rotations(n)] \
        + [p for _, _, p in affine(n)] + [hc[1]]
    rng2 = random.Random(seed * 77 + n)
    for _ in range(min(samples, 60)):
        p = list(range(n))
        rng2.shuffle(p)
        pool.append(tuple(p))
    rep["pool_size"] = len(pool)
    rep["H13I_violations_in_pool"] = sum(1 for p in pool if toric_I(p) > bound)
    rep["Pprime_violations_in_pool"] = sum(1 for p in pool if row_statistic(p) > n * bound)
    rep["minJ_violations_in_pool"] = sum(1 for p in pool if min_J(p) > bound)
    rep["max_I_minus_minJ_in_pool"] = str(max(Fraction(toric_I(p)) - min_J(p) for p in pool))
    rep["c_violations_in_pool"] = sum(
        1 for p in pool if c_invariant(p) > Fraction((n - 1) * (n - 2), 6))
    rep["c_equality_only_reflections"] = all(
        is_reflection(p)
        for p in pool if c_invariant(p) == Fraction((n - 1) * (n - 2), 6))
    return rep


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--ns", default="12,15,20,30")
    ap.add_argument("--samples", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--restarts", type=int, default=6)
    ap.add_argument("--sweeps", type=int, default=40)
    args = ap.parse_args(argv)
    out = {"tool": "toric_sampling-1.0", "seed": args.seed, "runs": []}
    for n in [int(x) for x in args.ns.split(",")]:
        rep = run(n, args.samples, args.seed, args.restarts, args.sweeps)
        out["runs"].append(rep)
        print(json.dumps(rep, ensure_ascii=False), flush=True)
    print(json.dumps({"summary": out["tool"], "ns": args.ns}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
