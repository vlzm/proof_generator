"""h13i_avg.py -- candidate route to H13-I: adaptive average over one cut coordinate.

H13-I (PLAN Sec 8, docs/notes/h13_line_model.md Sec 6): for every permutation pi of
Z_n, I(pi) = min_{p0,v0} inv(cut) <= floor((n-1)^2/4).

The uniform average of inv over all n^2 cuts (p0, v0) is already known to fail
(h13_line_model.md Sec 5: it is too high on e.g. sigma_n / identity-like inputs).
This script checks a DIFFERENT, not previously tried averaging (per the session's
starting notes): average only over p0, taking for each p0 the TRUE minimum over
v0, i.e.

    g(p0) = min_{v0} inv(p0, v0),          phi(pi) = sum_{p0=0}^{n-1} g(p0).

CANDIDATE LEMMA:  phi(pi) <= n * floor((n-1)^2/4)  for every pi.
If true for all n, H13-I follows by pigeonhole (some p0 has g(p0) <= average
<= floor((n-1)^2/4), and I(pi) <= g(p0)).

Two independent things are checked here, both purely combinatorial (no oracle
needed -- this is elementary permutation/inversion arithmetic, not an LRX-word
claim):

1. --verify-f: brute-force check of the known fact quoted in this session's
   task (f(k, m) = n(k+m) - 2km is the number of cuts (p0, v0) out of n^2 for
   which a pair with position gap k and value gap m is inverted).  Sanity
   check before building on it (rule: don't take given facts as given).

2. --search: NOT exhaustive.  Simulated-annealing / hill-climbing search for a
   permutation with phi(pi) > n * floor((n-1)^2/4), for n beyond exhaustive
   reach (12).  A search finding nothing is evidence, not a proof (rule 9,
   AGENTS.md); a search finding a violation is a genuine counterexample to the
   candidate lemma (not automatically to H13-I itself, since phi(pi) is only
   an upper-bound tool for I(pi)).

Exhaustive n <= 12 is done by the companion C program experiments/h13i_avg.c
(same recurrence, compiled and run separately; see data/runs/h13i_avg/).

Usage:
    python3 experiments/h13i_avg.py --verify-f --nmax 9
    python3 experiments/h13i_avg.py --search --n 20 50 100 --iters 4000 --restarts 6
Version h13i_avg-1.0 (session 9).
"""

import argparse
import itertools
import json
import math
import os
import random
import sys
import time

VERSION = "h13i_avg-1.0"
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "data", "runs", "h13i_avg")


def bound(n):
    return (n - 1) ** 2 // 4


def inv_count(seq):
    n = len(seq)
    c = 0
    for i in range(n):
        for j in range(i + 1, n):
            if seq[i] > seq[j]:
                c += 1
    return c


def verify_f(n, trials, rng):
    """Brute-force check: for random pairs (i, j) and random values (vi, vj),
    count cuts (p0, v0) in Z_n x Z_n for which the pair is inverted, and
    compare to f(k, m) = n(k+m) - 2km."""
    checked = 0
    for _ in range(trials):
        i = rng.randrange(n)
        j = rng.randrange(n)
        if i == j:
            continue
        k = (j - i) % n
        vi = rng.randrange(n)
        vj = rng.randrange(n)
        if vi == vj:
            continue
        m = (vj - vi) % n
        cnt = 0
        for p0 in range(n):
            for v0 in range(n):
                Pi = (i - p0) % n
                Pj = (j - p0) % n
                Vi = (vi - v0) % n
                Vj = (vj - v0) % n
                if (Pi < Pj and Vi > Vj) or (Pi > Pj and Vi < Vj):
                    cnt += 1
        f = n * (k + m) - 2 * k * m
        checked += 1
        if cnt != f:
            return False, checked, (i, j, k, m, cnt, f)
    return True, checked, None


def g_of_p0(line, n):
    """g(p0) via the telescoping recurrence inv(p0, v+1) - inv(p0, v)
    = n - 1 - 2*pos(v).  line = pi values at the n line positions for this p0."""
    inv0 = inv_count(line)
    pos = [0] * n
    for idx, val in enumerate(line):
        pos[val] = idx
    S = 0
    best = inv0
    for v in range(n - 1):
        D = n - 1 - 2 * pos[v]
        S += D
        cur = inv0 + S
        if cur < best:
            best = cur
    return best


def phi(pi, n):
    return sum(g_of_p0([pi[(p0 + j) % n] for j in range(n)], n) for p0 in range(n))


def hillclimb(n, iters, seed):
    rng = random.Random(seed)
    pi = list(rng.sample(range(n), n))
    cur = phi(pi, n)
    best = cur
    best_pi = pi[:]
    T0 = float(n)
    for it in range(iters):
        T = max(T0 * (1 - it / iters), 0.01)
        a, b = rng.sample(range(n), 2)
        pi[a], pi[b] = pi[b], pi[a]
        val = phi(pi, n)
        if val >= cur or rng.random() < math.exp((val - cur) / T):
            cur = val
        else:
            pi[a], pi[b] = pi[b], pi[a]
            continue
        if val > best:
            best = val
            best_pi = pi[:]
    return best, best_pi


def search_n(n, iters, restarts, seed0):
    reflection = [(-i) % n for i in range(n)]
    b_refl = phi(reflection, n)
    bnd = bound(n) * n
    result = {"n": n, "bound_times_n": bnd, "phi_reflection": b_refl,
              "reflection_status": "OK" if b_refl <= bnd else "VIOLATION",
              "restarts": []}
    overall_best = b_refl
    for r in range(restarts):
        t0 = time.time()
        b, p = hillclimb(n, iters, seed0 + r)
        dt = time.time() - t0
        result["restarts"].append({"seed": seed0 + r, "best": b, "time_s": round(dt, 2)})
        overall_best = max(overall_best, b)
    result["overall_best"] = overall_best
    result["status"] = "OK" if overall_best <= bnd else "VIOLATION"
    return result


def known_hard_families():
    """Known adversarial permutations from earlier H10/H12/H13 sessions
    (affine pi(i) = a*i + b, the exhaustively-found n = 10 H10 counterexample,
    reflections and single-transposition perturbations of reflections): check
    phi(pi) against n * floor((n-1)^2/4)."""
    cases = []
    cases.append(("H10 counterexample n=10", (9, 2, 5, 0, 1, 8, 7, 6, 3, 4), 10))
    cases.append(("affine n=17 (a=4,b=10)", [(4 * i + 10) % 17 for i in range(17)], 17))
    cases.append(("affine n=18 (a=5,b=16)", [(5 * i + 16) % 18 for i in range(18)], 18))
    cases.append(("affine n=54 (a=5,b=16)", [(5 * i + 16) % 54 for i in range(54)], 54))
    out = []
    for name, pi, n in cases:
        v = phi(list(pi), n)
        bnd = bound(n) * n
        out.append({"name": name, "n": n, "phi": v, "bound_times_n": bnd,
                     "status": "OK" if v <= bnd else "VIOLATION"})
    # affine sweep: all coprime (a, b) for n in {10, 18, 30}, b in {0, 1, n//2}
    sweep_viol = []
    for n in (10, 18, 30):
        for a in range(1, n):
            if math.gcd(a, n) != 1:
                continue
            for b in (0, 1, n // 2):
                pi = [(a * i + b) % n for i in range(n)]
                v = phi(pi, n)
                bnd = bound(n) * n
                if v > bnd:
                    sweep_viol.append({"n": n, "a": a, "b": b, "phi": v, "bound_times_n": bnd})
    # reflections and their single-transposition perturbations, n in {20, 30, 50}
    refl_out = []
    for n in (20, 30, 50):
        bnd = bound(n) * n
        refl = [(-i) % n for i in range(n)]
        base = phi(refl, n)
        worst = base
        for a in range(n):
            for b in range(a + 1, n):
                p = refl[:]
                p[a], p[b] = p[b], p[a]
                v = phi(p, n)
                if v > worst:
                    worst = v
        refl_out.append({"n": n, "phi_reflection": base, "bound_times_n": bnd,
                          "worst_single_transposition": worst,
                          "status": "OK" if worst <= bnd else "VIOLATION"})
    return {"cases": out, "affine_sweep_violations": sweep_viol, "reflection_perturbations": refl_out}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify-f", action="store_true")
    ap.add_argument("--nmax", type=int, default=9)
    ap.add_argument("--trials", type=int, default=300)
    ap.add_argument("--search", action="store_true")
    ap.add_argument("--families", action="store_true")
    ap.add_argument("--n", type=int, nargs="*", default=[20, 30, 50, 100])
    ap.add_argument("--iters", type=int, default=4000)
    ap.add_argument("--restarts", type=int, default=6)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    report = {"version": VERSION}

    if args.verify_f:
        rng = random.Random(args.seed)
        out = {}
        all_ok = True
        for n in range(4, args.nmax + 1):
            ok, checked, ce = verify_f(n, args.trials, rng)
            out[n] = {"ok": ok, "checked_pairs": checked, "counterexample": ce}
            all_ok = all_ok and ok
            print(f"verify_f n={n} checked={checked} ok={ok}" + ("" if ok else f" CE={ce}"))
        report["verify_f"] = out
        report["verify_f_all_ok"] = all_ok

    if args.search:
        out = []
        for n in args.n:
            r = search_n(n, args.iters, args.restarts, args.seed)
            out.append(r)
            print(f"search n={n} bound*n={r['bound_times_n']} phi(reflection)={r['phi_reflection']} "
                  f"overall_best={r['overall_best']} status={r['status']}")
        report["search"] = out

    if args.families:
        fam = known_hard_families()
        report["families"] = fam
        for c in fam["cases"]:
            print(f"family {c['name']}: phi={c['phi']} bound*n={c['bound_times_n']} status={c['status']}")
        print(f"affine sweep violations: {len(fam['affine_sweep_violations'])}")
        for r in fam["reflection_perturbations"]:
            print(f"reflection perturbations n={r['n']}: phi(refl)={r['phi_reflection']} "
                  f"worst_1swap={r['worst_single_transposition']} status={r['status']}")

    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(report, f, indent=2)
    print("wrote", os.path.join(OUT, "report.json"))


if __name__ == "__main__":
    main()
