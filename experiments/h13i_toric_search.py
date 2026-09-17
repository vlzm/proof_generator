"""H13-I (session 9): targeted local search for a counterexample at n >= 11.

Claim H13-I: I(pi) <= floor((n-1)^2/4) for every pi in S_n, where I(pi) is the
double-cut line-model quantity of experiments/line_model.py (min over q, c of
the inversion count of the relabelled line).  Exhaustive brute force over all
n! permutations is only cheap enough up to n = 12 (experiments/line_toric_max_i.c,
O(n^3) per permutation); this script instead runs simulated annealing that
tries to MAXIMIZE I(pi), starting from random permutations and from perturbed
reflections, for n where exhaustive search is not affordable (n >= 12; also
used to cross-check n = 11..12 against the exhaustive result).

I(pi) is computed directly from the definition (O(n^2) cuts * O(n^2)
inversions = O(n^4) per permutation); at n = 15 this is ~1.4 ms/call
(measured), so a few hundred thousand evaluations per n are affordable in
minutes.  This does NOT reuse experiments/line_toric_max_i.c or
experiments/line_model.py: it is a fresh, independent evaluator, used only
for local search, not as a certificate (the exhaustive C program and the
distance tables remain the ground truth for n <= 11/12).

Usage: python3 experiments/h13i_toric_search.py N ITERS TRIALS [--seed S]
Output: best I(pi) found and an argmax; prints one line per trial.
Version h13i_toric_search-1.0.
"""

import argparse
import random
import time


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def I_of(pi):
    n = len(pi)
    best = None
    for q in range(n):
        for c in range(n):
            shift = (q + 1 - c) % n
            w = [(pi[(q + 1 + j) % n] - shift) % n for j in range(n)]
            inv = inversions(w)
            if best is None or inv < best:
                best = inv
    return best


def reflection(n, h):
    return tuple((h - i) % n for i in range(n))


def anneal(n, iters, seed, start=None, t0=3.0, t1=0.02):
    random.seed(seed)
    pi = list(range(n))
    random.shuffle(pi)
    if start is not None:
        pi = list(start)
    cur = I_of(pi)
    best, best_pi = cur, list(pi)
    for it in range(iters):
        temp = t0 * (t1 / t0) ** (it / iters)
        i, j = random.sample(range(n), 2)
        pi[i], pi[j] = pi[j], pi[i]
        val = I_of(pi)
        if val >= cur or random.random() < pow(2.718281828, (val - cur) / temp):
            cur = val
            if val > best:
                best, best_pi = val, list(pi)
        else:
            pi[i], pi[j] = pi[j], pi[i]
    return best, best_pi


def run(n, iters, trials, seed0):
    target = (n - 1) ** 2 // 4
    print(f"n={n} target floor((n-1)^2/4)={target} iters={iters} trials={trials}")
    overall_best, overall_pi = -1, None
    t0 = time.time()
    for tr in range(trials):
        random.seed(seed0 + tr)
        if tr % 2 == 0:
            start = None
        else:
            h = random.randrange(n)
            start = list(reflection(n, h))
            for _ in range(random.randint(0, 3)):
                i, j = random.sample(range(n), 2)
                start[i], start[j] = start[j], start[i]
        best, best_pi = anneal(n, iters, seed=seed0 + tr, start=start)
        print(f"  trial {tr}: best I = {best}  pi={best_pi}")
        if best > overall_best:
            overall_best, overall_pi = best, best_pi
    dt = time.time() - t0
    print(f"n={n}: OVERALL BEST I = {overall_best} (target {target})  time={dt:.1f}s")
    print("argmax:", overall_pi)
    return {"n": n, "iters": iters, "trials": trials, "seed0": seed0, "target": target,
            "best": overall_best, "argmax": overall_pi, "exceeds_target": overall_best > target,
            "seconds": round(dt, 1)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("iters", type=int)
    ap.add_argument("trials", type=int)
    ap.add_argument("--seed", type=int, default=1000)
    args = ap.parse_args()
    run(args.n, args.iters, args.trials, args.seed)


if __name__ == "__main__":
    main()
