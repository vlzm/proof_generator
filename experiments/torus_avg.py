"""H13-I candidate route: average of the per-rotation optimum (session 9).

I(pi) = min_{a,b in Z_n} inv(w_{a,b}) is the minimum, over all n^2 double cuts
(position rotation a, value rotation b), of the number of inversions of the
relabelled line (docs/notes/h13_line_model.md).  Averaging inv(w_{a,b})
uniformly over all n^2 cuts does NOT bound I(pi) by floor((n-1)^2/4) (PLAN
S8, "среднее ... ~n^2/3 для id"): the naive Markov argument is too weak.

This script tests a two-step route instead.  Fix a; define
    M(a) = min_b inv(w_{a,b})
(the best value-rotation for that one position-rotation).  Trivially
I(pi) = min_a M(a) <= mean_a M(a), so it is enough to bound the MEAN of M(a)
over the n position-rotations, not the min directly.

Closed form for M(a) (telescoping recursion, derived this session).  Fix the
position-rotation a and write x = (pi(a), pi(a+1 mod n), ..., pi(a+n-1 mod n))
for the resulting line.  As b increases by 1, exactly one element (the one
holding shifted-value 0 under b, i.e. original value b, at line-position p)
flips from "smallest" to "largest" under the b -> b+1 relabelling, and every
pair involving it flips inverted/not-inverted; no other pair's order changes.
Hence
    inv(x, b+1) = inv(x, b) + (n - 1 - 2 p_b),   p_b = position of value b in x.
Summing from b=0 (writing rho = x^{-1}, P(k) = rho(0)+...+rho(k-1)):
    inv(x, b) = inv(x, 0) + b*(n-1) - 2*P(b),   0 <= b <= n  (P(n) = n(n-1)/2
    closes the cycle: inv(x, n) = inv(x, 0)).
So M(a) = min_{0<=k<=n} [inv(x,0) + k*(n-1) - 2*P(k)], computable in O(n) after
O(n^2) (or O(n log n) with a Fenwick tree; not needed at these n) to get
inv(x,0) and the p_b array.  Cross-checked against brute force min_b (all b,
naive O(n) x O(n^2)) and against the full O(n^2) two-cut search for I(pi) on
random inputs and on the two adversarial lines from h13_line_model.md S4/S1.7
(the exhaustive test below also re-derives I(pi) this way and can be diffed
against experiments/line_profile.c's I_n{n}.bin).

Conjecture tested here (H13-I''): mean_a M(a) <= floor((n-1)^2/4) for every
pi in S_n, with equality iff pi is a reflection (pi(i) = h - i mod n for some
h).  This is STRONGER than what H13-I needs (it implies I(pi) <= floor((n-1)^2
/4) since min <= mean) and appears more tractable: it reduces to bounding a
sum of n closed-form expressions instead of a joint minimum over n^2 points.
Not proved; this script only gathers exhaustive/sampled evidence.

Usage:
  python3 experiments/torus_avg.py --nmax 9            # exhaustive 2 <= n <= 9
  python3 experiments/torus_avg.py --sample 20,30,50 --trials 4000
Output: data/runs/torus_avg/report.json, report.md.  Version torus_avg-1.0.
"""

import argparse
import itertools
import json
import os
import random
import time

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "runs", "torus_avg")
VERSION = "torus_avg-1.0"


def inv_count(seq):
    n = len(seq)
    return sum(1 for i in range(n) for j in range(i + 1, n) if seq[i] > seq[j])


def M_fast(x):
    """min_b inv(x rotated in value by b), via the telescoping closed form."""
    n = len(x)
    inv0 = inv_count(x)
    p = [0] * n
    for pos, v in enumerate(x):
        p[v] = pos
    P = 0
    best = inv0  # k = 0 term
    for k in range(1, n):
        P += p[k - 1]
        g = inv0 + k * (n - 1) - 2 * P
        if g < best:
            best = g
    return best


def M_brute(x):
    n = len(x)
    return min(inv_count([(v - b) % n for v in x]) for b in range(n))


def rotations(pi):
    n = len(pi)
    return [tuple(pi[(a + i) % n] for i in range(n)) for a in range(n)]


def sum_M(pi):
    return sum(M_fast(x) for x in rotations(pi))


def I_of(pi):
    return min(M_fast(x) for x in rotations(pi))


def self_check(trials=3000, seed=1):
    """M_fast vs M_brute vs full two-cut search, small random inputs."""
    rng = random.Random(seed)
    for _ in range(trials):
        n = rng.randint(2, 11)
        x = list(range(n))
        rng.shuffle(x)
        a, b = M_brute(x), M_fast(x)
        if a != b:
            raise AssertionError(f"M_fast mismatch n={n} x={x}: brute={a} fast={b}")
    for _ in range(200):
        n = rng.randint(2, 8)
        pi = list(range(n))
        rng.shuffle(pi)
        fast_I = I_of(tuple(pi))
        brute_I = min(
            inv_count([(v - b) % n for v in [pi[(a + i) % n] for i in range(n)]])
            for a in range(n)
            for b in range(n)
        )
        if fast_I != brute_I:
            raise AssertionError(f"I mismatch n={n} pi={pi}: brute={brute_I} fast={fast_I}")
    return True


def target(n):
    return (n - 1) ** 2 // 4


def exhaustive(nmin, nmax):
    rows = []
    for n in range(nmin, nmax + 1):
        t0 = time.time()
        tgt = target(n)
        worst_avg = -1.0
        worst_avg_x = None
        ties_avg = 0
        count = 0
        for perm in itertools.permutations(range(n)):
            s = sum_M(perm)
            avg = s / n
            count += 1
            if avg > worst_avg + 1e-12:
                worst_avg = avg
                worst_avg_x = perm
                ties_avg = 1
            elif abs(avg - worst_avg) < 1e-9:
                ties_avg += 1
        elapsed = time.time() - t0
        ok = worst_avg <= tgt + 1e-9
        rows.append(
            {
                "n": n,
                "target": tgt,
                "max_mean_M": worst_avg,
                "argmax": list(worst_avg_x),
                "ties_at_max": ties_avg,
                "n_factorial": count,
                "ok": ok,
                "seconds": elapsed,
            }
        )
        print(f"n={n} target={tgt} max_mean_M={worst_avg} ties={ties_avg} "
              f"{'OK' if ok else 'VIOLATION'} ({elapsed:.1f}s, {count} perms)")
    return rows


def sampled(ns, trials, seed=12345):
    rng = random.Random(seed)
    rows = []
    for n in ns:
        tgt = target(n)
        candidates = []
        perturb_positions = range(n - 1) if n <= 20 else rng.sample(range(n - 1), 20)
        for h in range(n):
            refl = [(h - i) % n for i in range(n)]
            candidates.append(tuple(refl))
            for i in perturb_positions:
                r2 = refl[:]
                r2[i], r2[i + 1] = r2[i + 1], r2[i]
                candidates.append(tuple(r2))
        for _ in range(trials):
            p = list(range(n))
            rng.shuffle(p)
            candidates.append(tuple(p))
        worst_avg = -1.0
        worst_x = None
        for perm in candidates:
            avg = sum_M(perm) / n
            if avg > worst_avg:
                worst_avg = avg
                worst_x = perm
        ok = worst_avg <= tgt + 1e-9
        rows.append(
            {
                "n": n,
                "target": tgt,
                "max_mean_M_sampled": worst_avg,
                "n_candidates": len(candidates),
                "ok": ok,
            }
        )
        print(f"n={n} target={tgt} max_mean_M(sampled+reflections)={worst_avg} "
              f"{'OK' if ok else 'VIOLATION'} ({len(candidates)} candidates)")
    return rows


def check_single_b(nmin, nmax):
    """REFUTED (session 9, §7.2): is M(x) = min_b inv(x,b) <= floor((n-1)^2/4)
    for EVERY string x (fixed position order, i.e. a = 0 only, no rotation of
    positions)?  Exhaustive over all n! strings, not just rotations of a
    cyclic pi."""
    rows = []
    for n in range(nmin, nmax + 1):
        tgt = target(n)
        worst = -1
        worst_x = None
        for perm in itertools.permutations(range(n)):
            m = M_fast(list(perm))
            if m > worst:
                worst = m
                worst_x = perm
        ok = worst <= tgt
        rows.append({"n": n, "target": tgt, "max_M_a0_only": worst, "argmax": list(worst_x), "ok": ok})
        print(f"n={n} target={tgt} max M(x), a=0 only = {worst} argmax={worst_x} "
              f"{'OK' if ok else 'VIOLATION (refutes single-rotation simplification)'}")
    return rows


def check_diagonal(nmin, nmax):
    """REFUTED (session 9, §7.3): is min over the diagonal family
    (a, b=pi(a)) (n candidates instead of n^2) <= floor((n-1)^2/4)?"""
    rows = []
    for n in range(nmin, nmax + 1):
        tgt = target(n)
        worst = -1
        worst_x = None
        for perm in itertools.permutations(range(n)):
            best = None
            for a in range(n):
                x = [(perm[(a + j) % n] - perm[a]) % n for j in range(n)]
                c = inv_count(x)
                if best is None or c < best:
                    best = c
            if best > worst:
                worst = best
                worst_x = perm
        ok = worst <= tgt
        rows.append({"n": n, "target": tgt, "max_diagonal_min": worst, "argmax": list(worst_x), "ok": ok})
        print(f"n={n} target={tgt} max diagonal min = {worst} argmax={worst_x} "
              f"{'OK' if ok else 'VIOLATION (refutes diagonal-family simplification)'}")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=2)
    ap.add_argument("--nmax", type=int, default=9)
    ap.add_argument("--sample", type=str, default="")
    ap.add_argument("--trials", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--mode", choices=["mean_M", "single_b", "diagonal"], default="mean_M",
                     help="mean_M: H13-I'' (this script's main route); "
                          "single_b / diagonal: reproduce the two refuted simplifications (S7.2, S7.3)")
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    print("self-check (M_fast vs brute, I via M_fast vs full two-cut search)...")
    self_check()
    print("self-check passed.")

    if args.mode == "single_b":
        report = {"version": VERSION, "self_check": "passed", "single_b_refuted": check_single_b(args.nmin, args.nmax)}
        with open(os.path.join(OUT, "report.json"), "a") as f:
            f.write(json.dumps(report) + "\n")
        return
    if args.mode == "diagonal":
        report = {"version": VERSION, "self_check": "passed", "diagonal_refuted": check_diagonal(args.nmin, args.nmax)}
        with open(os.path.join(OUT, "report.json"), "a") as f:
            f.write(json.dumps(report) + "\n")
        return

    report = {"version": VERSION, "self_check": "passed"}
    if args.sample:
        ns = [int(x) for x in args.sample.split(",") if x.strip()]
        report["sampled"] = sampled(ns, args.trials, args.seed)
    else:
        report["exhaustive"] = exhaustive(args.nmin, args.nmax)

    with open(os.path.join(OUT, "report.json"), "a") as f:
        f.write(json.dumps(report) + "\n")
    print("appended to", os.path.join(OUT, "report.json"))


if __name__ == "__main__":
    main()
