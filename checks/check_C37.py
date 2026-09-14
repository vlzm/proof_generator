#!/usr/bin/env python3
"""Check C37: exact recursion identity for value-shift inversions, and
counterexamples showing that a single free shift (or two independent
single-shift strategies) is not enough to prove H13-I in general.

Usage:
    python3 checks/check_C37.py --identity-max 8
    python3 checks/check_C37.py --strategy-max 9
    python3 checks/check_C37.py --identity-max 8 --strategy-max 9 --sample-n 20,50,100
"""
import argparse
import itertools
import json
import os
import random
import time


def inv_count(seq):
    n = len(seq)
    c = 0
    for i in range(n):
        vi = seq[i]
        for j in range(i + 1, n):
            if vi > seq[j]:
                c += 1
    return c


def cross_inv(tau, c, pos):
    # X_c(tau) = #{(x,y): x < c <= y, pos(y) < pos(x)}
    cnt = 0
    for x in range(c):
        px = pos[x]
        for y in range(c, len(tau)):
            if pos[y] < px:
                cnt += 1
    return cnt


def shift(tau, c):
    n = len(tau)
    return [(x - c) % n for x in tau]


def positions(tau):
    n = len(tau)
    pos = [0] * n
    for idx, v in enumerate(tau):
        pos[v] = idx
    return pos


def check_identity(nmax, sample_ns, seed=0):
    report = {"exhaustive": {}, "sampled": {}}
    for n in range(4, nmax + 1):
        base_inv_cache = {}
        checked = 0
        for tau in itertools.permutations(range(n)):
            tau = list(tau)
            pos = positions(tau)
            base = inv_count(tau)
            for c in range(n):
                lhs = inv_count(shift(tau, c))
                rhs = base + c * (n - c) - 2 * cross_inv(tau, c, pos)
                assert lhs == rhs, (n, tau, c, lhs, rhs)
                checked += 1
        report["exhaustive"][n] = {"permutations": n and __import__("math").factorial(n), "pairs_checked": checked}
        print(f"identity exhaustive n={n}: OK ({checked} (tau,c) pairs)")

    rng = random.Random(seed)
    for n in sample_ns:
        trials = 500
        for _ in range(trials):
            tau = list(range(n))
            rng.shuffle(tau)
            pos = positions(tau)
            base = inv_count(tau)
            for c in range(n):
                lhs = inv_count(shift(tau, c))
                rhs = base + c * (n - c) - 2 * cross_inv(tau, c, pos)
                assert lhs == rhs, (n, tau, c, lhs, rhs)
        report["sampled"][n] = {"trials": trials}
        print(f"identity sampled n={n}: OK ({trials} random permutations x {n} shifts)")
    return report


def best_c(tau):
    n = len(tau)
    return min(inv_count(shift(tau, c)) for c in range(n))


def inverse(pi):
    n = len(pi)
    inv = [0] * n
    for i, v in enumerate(pi):
        inv[v] = i
    return inv


def check_strategy_failures(nmax):
    report = {}
    for n in range(4, nmax + 1):
        bound = ((n - 1) ** 2) // 4
        worst_q0 = 0
        worst_combo = 0
        fail_q0_example = None
        fail_combo_example = None
        t0 = time.time()
        for pi in itertools.permutations(range(n)):
            pi = list(pi)
            a = best_c(pi)
            if a > worst_q0:
                worst_q0 = a
            if a > bound and fail_q0_example is None:
                fail_q0_example = pi
            b = best_c(inverse(pi))
            combo = min(a, b)
            if combo > worst_combo:
                worst_combo = combo
            if combo > bound and fail_combo_example is None:
                fail_combo_example = pi
        dt = time.time() - t0
        report[n] = {
            "bound_floor_n_1_sq_4": bound,
            "max_q0_only": worst_q0,
            "max_combo_q0_or_c0": worst_combo,
            "q0_only_fails": worst_q0 > bound,
            "combo_fails": worst_combo > bound,
            "min_n_witness_q0_only": list(fail_q0_example) if fail_q0_example else None,
            "min_n_witness_combo": list(fail_combo_example) if fail_combo_example else None,
            "seconds": round(dt, 1),
        }
        print(f"n={n}: bound={bound} max(q=0 only)={worst_q0} max(combo A,B)={worst_combo} "
              f"[{dt:.1f}s]")
    return report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--identity-max", type=int, default=8)
    ap.add_argument("--strategy-max", type=int, default=9)
    ap.add_argument("--sample-n", type=str, default="20,50")
    ap.add_argument("--out", type=str, default="data/runs/check_C37/report.json")
    args = ap.parse_args()

    sample_ns = [int(x) for x in args.sample_n.split(",") if x.strip()]

    identity_report = check_identity(args.identity_max, sample_ns)
    strategy_report = check_strategy_failures(args.strategy_max)

    out = {"identity": identity_report, "strategy_failures": strategy_report}
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print("PASS")
    print(f"report written to {args.out}")


if __name__ == "__main__":
    main()
