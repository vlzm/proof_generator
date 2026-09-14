"""h13i_shift_candidates.py -- exploring candidate restrictions for H13-I.

H13-I (PLAN §8, docs/notes/h13_line_model.md §6): every permutation pi of
Z_n has a double cut (position rotation q, value rotation c) with
inv(w) <= floor((n-1)^2/4), where w_j = pi(q+1+j) - (q+1-c) mod n.

This script tests whether the position rotation can be restricted to a
small, n-independent (or slowly growing) candidate set without losing the
bound floor((n-1)^2/4), i.e. whether

    min_{r in R} min_{s in Z_n} inv(shift(pi, r, s)) <= floor((n-1)^2/4)

holds for all pi, for a fixed small set R (e.g. R = {0}, R = {0, 1}) or for
R = {(i0, pi(i0)) : i0} (cut through a fixed point of pi).

Modes:
  exhaustive-r0     all pi, n in [nmin, nmax]: r fixed at 0, best s only.
  exhaustive-diag   all pi, n in [nmin, nmax]: r = i0, s = pi(i0), best i0.
  exhaustive-rset   all pi, n in [nmin, nmax]: r in RSET (--rset "0,1"), best s.
  local-search-rset hill-climbing over random pi, n in [nmin, nmax]: same
                    restricted score as exhaustive-rset, to look for growing
                    counterexamples beyond the exhaustive range.

Version h13i_shift_candidates-1.0.
"""
import argparse
import itertools
import json
import random
import time


def inv_count(seq):
    n = len(seq)
    c = 0
    for j in range(n):
        vj = seq[j]
        for k in range(j + 1, n):
            if vj > seq[k]:
                c += 1
    return c


def shifted(pi, n, r, s):
    return [(pi[(r + j) % n] - s) % n for j in range(n)]


def score_rset(pi, n, rset):
    best = None
    for r in rset:
        for s in range(n):
            m = inv_count(shifted(pi, n, r, s))
            if best is None or m < best:
                best = m
    return best


def score_diag(pi, n):
    best = None
    for i0 in range(n):
        m = inv_count(shifted(pi, n, i0, pi[i0]))
        if best is None or m < best:
            best = m
    return best


def bound(n):
    return ((n - 1) ** 2) // 4


def run_exhaustive(nmin, nmax, mode, rset):
    for n in range(nmin, nmax + 1):
        b = bound(n)
        worst = -1
        worst_pi = None
        exceed = 0
        checked = 0
        t0 = time.time()
        for pi in itertools.permutations(range(n)):
            checked += 1
            if mode == "exhaustive-r0":
                m = score_rset(pi, n, [0])
            elif mode == "exhaustive-rset":
                m = score_rset(pi, n, rset)
            elif mode == "exhaustive-diag":
                m = score_diag(pi, n)
            else:
                raise ValueError(mode)
            if m > b:
                exceed += 1
            if m > worst:
                worst = m
                worst_pi = pi
        dt = time.time() - t0
        print(json.dumps({
            "mode": mode, "n": n, "rset": rset if mode == "exhaustive-rset" else None,
            "bound": b, "worst": worst, "excess": worst - b,
            "exceed_count": exceed, "checked": checked,
            "worst_pi": list(worst_pi), "seconds": round(dt, 2),
        }))


def run_local_search(nmin, nmax, rset, iters, restarts, seed):
    rng = random.Random(seed)
    for n in range(nmin, nmax + 1):
        b = bound(n)
        global_best = -1
        global_pi = None
        for _ in range(restarts):
            pi = list(range(n))
            rng.shuffle(pi)
            cur = score_rset(pi, n, rset)
            for _ in range(iters):
                i, j = rng.sample(range(n), 2)
                pi[i], pi[j] = pi[j], pi[i]
                new = score_rset(pi, n, rset)
                if new >= cur:
                    cur = new
                else:
                    pi[i], pi[j] = pi[j], pi[i]
            if cur > global_best:
                global_best = cur
                global_pi = pi[:]
        print(json.dumps({
            "mode": "local-search-rset", "n": n, "rset": rset,
            "bound": b, "best_found": global_best, "excess": global_best - b,
            "iters": iters, "restarts": restarts, "seed": seed,
            "best_pi": global_pi,
            "note": "heuristic lower bound on the true worst case, not exhaustive",
        }))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=[
        "exhaustive-r0", "exhaustive-diag", "exhaustive-rset", "local-search-rset",
    ])
    ap.add_argument("--nmin", type=int, default=4)
    ap.add_argument("--nmax", type=int, default=9)
    ap.add_argument("--rset", type=str, default="0,1")
    ap.add_argument("--iters", type=int, default=2500)
    ap.add_argument("--restarts", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rset = [int(x) for x in args.rset.split(",") if x != ""]

    if args.mode == "local-search-rset":
        run_local_search(args.nmin, args.nmax, rset, args.iters, args.restarts, args.seed)
    else:
        run_exhaustive(args.nmin, args.nmax, args.mode, rset)


if __name__ == "__main__":
    main()
