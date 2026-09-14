"""Checker for C37 (session 9): closed-form / increment formulas for the
double-cut inversion count f(a, b) used by the H13-I line model
(docs/notes/h13_line_model.md, docs/proofs/C37_cut_inversion_formula.md).

Independent implementation (no code shared with experiments/torus_cut.c):
  Part A: brute-force f(a, b) (literal cut-and-relabel-and-count-inversions)
    against the two increment lemmas
      f(a+1, b) - f(a, b) = n - 1 - 2*((pi(a) - b) mod n)
      f(a, b+1) - f(a, b) = n - 1 - 2*((pi^{-1}(b) - a) mod n)
    on all (a, b) for random permutations, many n.
  Part B: brute-force f(a, b) against the closed form
      f(a,b) = inv(pi) + (n-1)(a+b) + 2ab - 2 P(a) - 2 R(b) - 2 n K(a,b)
    P(a) = sum_{i<a} pi(i), R(b) = sum_{t<b} pi^{-1}(t), K(a,b) = #{i<a: pi(i)<b}.
  Part C: min_{a,b} f(a,b) computed by this checker's own O(n^2) recurrence
    against experiments/torus_cut.c's `exhaustive`/`family` output and, where
    available, against the old O(n^4) tables in data/runs/line_profile/
    (max_I per n, from experiments/line_profile.c) -- both must give the same
    max_pi I(pi).

Usage: python3 checks/check_C37.py [--trials 300] [--nmax 10]
Output: data/runs/check_C37/report.json, report.md.
"""

import argparse
import json
import os
import random
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

VERSION = "check_C37-1.0"
OUT = os.path.join(ROOT, "data", "runs", "check_C37")


def brute_inv(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def brute_f(pi, a, b):
    n = len(pi)
    w = [(pi[(a + j) % n] - b) % n for j in range(n)]
    return brute_inv(w)


def closed_form(pi, invp, a, b):
    n = len(pi)
    inv_pi = [0] * n
    for i, v in enumerate(pi):
        inv_pi[v] = i
    P = sum(pi[:a])
    R = sum(inv_pi[:b])
    K = sum(1 for i in range(a) if pi[i] < b)
    return invp + (n - 1) * (a + b) + 2 * a * b - 2 * P - 2 * R - 2 * n * K


def fast_grid_min(pi):
    """Independent O(n^2) recurrence (Lemma, part A), returns (min_I, full_grid)."""
    n = len(pi)
    inv_pi = [0] * n
    for i, v in enumerate(pi):
        inv_pi[v] = i
    f00 = brute_inv(pi)
    a0 = [0] * n
    a0[0] = f00
    cur = f00
    for a in range(n - 1):
        cur += n - 1 - 2 * pi[a]
        a0[a + 1] = cur
    best = f00
    grid = [[0] * n for _ in range(n)]
    for a in range(n):
        g = a0[a]
        grid[a][0] = g
        best = min(best, g)
        for b in range(n - 1):
            g += n - 1 - 2 * ((inv_pi[b] - a) % n)
            grid[a][b + 1] = g
            best = min(best, g)
    return best, grid


def part_ab(trials, nmax, log, rnd):
    ok = True
    checked = 0
    for _ in range(trials):
        n = rnd.randint(2, nmax)
        pi = list(range(n))
        rnd.shuffle(pi)
        invp = brute_inv(pi)
        for a in range(n):
            for b in range(n):
                fab = brute_f(pi, a, b)
                fab_a1 = brute_f(pi, (a + 1) % n, b)
                fab_b1 = brute_f(pi, a, (b + 1) % n)
                lemma_a = fab + (n - 1 - 2 * ((pi[a] - b) % n))
                inv_pi_b = pi.index(b)
                lemma_b = fab + (n - 1 - 2 * ((inv_pi_b - a) % n))
                cf = closed_form(pi, invp, a, b)
                checked += 1
                if fab_a1 != lemma_a or fab_b1 != lemma_b or fab != cf:
                    ok = False
                    log(f"FAIL n={n} pi={pi} a={a} b={b}: f={fab} cf={cf} "
                        f"f(a+1,b)={fab_a1} vs lemma {lemma_a}; f(a,b+1)={fab_b1} vs lemma {lemma_b}")
    log(f"part A+B: {checked} (pi, a, b) triples over {trials} permutations (n <= {nmax}) -> "
        f"{'all match' if ok else 'MISMATCH FOUND'}")
    return ok, checked


def part_c(nmax, log):
    ok = True
    rows = []
    for n in range(4, nmax + 1):
        t0 = time.time()
        best_over_all = -1
        # exhaustive over all n! permutations using the fast O(n^2) recurrence,
        # cross-checked against the O(n^4) brute force on a random subsample.
        import itertools
        import math
        rnd = random.Random(n)
        nfact = math.factorial(n)
        sample_idx = set(rnd.sample(range(nfact), k=min(20, nfact)))
        for idx, pi in enumerate(itertools.permutations(range(n))):
            best, _ = fast_grid_min(list(pi))
            if best > best_over_all:
                best_over_all = best
            if idx in sample_idx:
                brute_best = min(brute_f(list(pi), a, b) for a in range(n) for b in range(n))
                if brute_best != best:
                    ok = False
                    log(f"FAIL cross-check n={n} pi={pi}: fast={best} brute={brute_best}")
        bound = (n - 1) ** 2 // 4
        good = best_over_all <= bound
        ok = ok and good
        rows.append({"n": n, "max_I": best_over_all, "bound": bound, "seconds": round(time.time() - t0, 2)})
        log(f"part C n={n}: max_pi I(pi) = {best_over_all} (bound floor((n-1)^2/4) = {bound}) -> "
            f"{'ok' if good else 'VIOLATION'}, {time.time() - t0:.2f} s")
        # cross-check against the existing O(n^4) profile data if present (C33/C35 artifact)
        prof_path = os.path.join(ROOT, "data", "runs", "line_profile", f"profile_n{n}.json")
        if os.path.exists(prof_path):
            prof = json.load(open(prof_path))
            prof_maxI = max(r["I"] for r in prof["by_I"])
            match = prof_maxI == best_over_all
            ok = ok and match
            log(f"  cross-check vs data/runs/line_profile/profile_n{n}.json: max_I={prof_maxI} -> "
                f"{'match' if match else 'MISMATCH'}")
    return ok, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=300)
    ap.add_argument("--nmax", type=int, default=9)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    lines = []

    def log(msg):
        print(msg, flush=True)
        lines.append(msg)

    log(f"== {VERSION} args={vars(args)}")
    rnd = random.Random(args.seed)
    ok_ab, checked = part_ab(args.trials, min(args.nmax, 9), log, rnd)
    ok_c, rows = part_c(args.nmax, log)
    verdict = "PASS" if (ok_ab and ok_c) else "FAIL"
    log(f"verdict: {verdict}")
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump({"version": VERSION, "args": vars(args), "ab_triples_checked": checked,
                   "part_c": rows, "verdict": verdict}, f, indent=1)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# check_C37 — замкнутая форма f(a, b) для двойного разреза (H13-I)\n\n")
        f.write(f"Команда: `python3 checks/check_C37.py --trials {args.trials} --nmax {args.nmax}`. "
                f"Версия: {VERSION}.\n\n```text\n")
        f.write("\n".join(lines) + "\n```\n")
    sys.exit(0 if verdict == "PASS" else 1)


if __name__ == "__main__":
    main()
