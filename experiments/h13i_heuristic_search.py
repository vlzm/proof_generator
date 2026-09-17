"""h13i_heuristic_search.py — simulated-annealing search for a counterexample
to H13-I: I(pi) <= floor((n-1)^2/4) for every permutation pi of {0,...,n-1},
where I(pi) = min over double toric cuts (a, b) of the linear inversion count
of the relabelled line (see docs/notes/h13_line_model.md).

Not exhaustive: this is a search for errors (AGENTS.md rule 9), not a proof.
Use experiments/h13i_exhaustive.c for exhaustive verification (feasible up to
n ~ 12-13 thanks to the O(n^2)-per-permutation incremental formula below).

I(pi) is computed via the double-incremental recurrence:
  - rotating positions by 1 (move front element to back) changes inv by
    (n - 1 - 2*w_0), w_0 = value currently at the front;
  - rotating values by 1 (subtract 1 mod n from every value) changes inv by
    (n - 1 - 2*j*), j* = position of the value about to become 0.
Verified against brute force over all n^2 cuts for n = 3..8 (exact match).

Version h13i_heuristic_search-1.0.
Usage: python3 h13i_heuristic_search.py [n1 n2 ...]
"""
import random
import math
import sys
import time
import json


def I_fast(pi):
    n = len(pi)
    w = pi[:]
    inv00 = sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])
    best = inv00
    inv_a = inv00
    rot = pi[:]
    for a in range(n):
        if a > 0:
            front = rot[0]
            rot = rot[1:] + [front]
            inv_a += n - 1 - 2 * front
        if inv_a < best:
            best = inv_a
        pos_of_value = [0] * n
        for j, v in enumerate(rot):
            pos_of_value[v] = j
        inv_b = inv_a
        for b in range(n - 1):
            jstar = pos_of_value[b]
            inv_b += n - 1 - 2 * jstar
            if inv_b < best:
                best = inv_b
    return best


def target(n):
    return ((n - 1) * (n - 1)) // 4


def reflection(n, h=0):
    return [(h - i) % n for i in range(n)]


def anneal(n, iters, restarts, seed, reflection_seeds=True):
    rng = random.Random(seed)
    global_best = -1
    global_pi = None
    for r in range(restarts):
        if reflection_seeds and r < min(n, 4):
            pi = reflection(n, r)
        else:
            pi = list(range(n))
            rng.shuffle(pi)
        cur = I_fast(pi)
        best_local = cur
        t0 = 3.0
        for it in range(iters):
            temp = t0 * (1 - it / iters) + 0.05
            i, j = rng.sample(range(n), 2)
            pi[i], pi[j] = pi[j], pi[i]
            new = I_fast(pi)
            delta = new - cur
            if delta >= 0 or rng.random() < math.exp(delta / temp):
                cur = new
                if cur > best_local:
                    best_local = cur
            else:
                pi[i], pi[j] = pi[j], pi[i]
        if best_local > global_best:
            global_best = best_local
            global_pi = pi[:]
    return global_best, global_pi


def main():
    ns = [int(x) for x in sys.argv[1:]] or [11, 12, 13, 14, 16, 18, 20, 24, 28, 32]
    results = []
    for n in ns:
        t0 = time.time()
        tgt = target(n)
        best, pi = anneal(n, iters=1500 if n < 20 else 800, restarts=5, seed=0)
        dt = time.time() - t0
        holds = best <= tgt
        results.append({
            "n": n, "target": tgt, "best_found_I": best,
            "H13_I_violated": not holds, "argmax": pi, "seconds": dt,
        })
        print(json.dumps(results[-1]))
    return results


if __name__ == "__main__":
    main()
