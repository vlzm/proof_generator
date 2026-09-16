"""H13-I (toric-class inversions): does a *single* free cut parameter suffice?

Context: I(pi) = min over all n^2 double cuts (a, b) of inv(w_{a,b}), where
w_{a,b}(j) = (pi((a+j) mod n) - b) mod n, j = 0..n-1 (position rotation a,
value rotation b).  H13-I conjectures max_pi I(pi) = floor((n-1)^2/4)
(VERIFIED 4 <= n <= 10 by experiments/line_profile.c).  Open: a proof.

This script tests whether the full 2-parameter search over (a, b) can be
replaced by a search over only ~n candidates (one free parameter), which
would make a direct extremal/inductive proof much easier.  Two families:

1. Fix the position cut a = 0 and vary only the value cut b (n candidates).
2. Couple b to a by an explicit rule b = f(pi, a) (n candidates), for three
   natural rules that make the line start at a fixed value.

For each family we report, over pi, the worst case of
    min_{candidate} inv(w)
compared against the bound floor((n-1)^2/4) and against the true I(pi)
(computed by brute force over all n^2 cuts, feasible up to n ~ 9).

Usage:
    python3 experiments/h13i_coupling.py --nmax 8
    python3 experiments/h13i_coupling.py --n 9 --sample 20000 --seed 1
    python3 experiments/h13i_coupling.py --n 20 --sample 5000 --seed 1

Output: data/runs/h13i_coupling/report.md (appended per run block).
Version h13i_coupling-1.0 (session 9).
"""

import argparse
import itertools
import json
import os
import random
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def inv_count(seq):
    n = len(seq)
    c = 0
    for i in range(n):
        si = seq[i]
        for j in range(i + 1, n):
            if si > seq[j]:
                c += 1
    return c


def rotated(perm, a):
    n = len(perm)
    return perm[a:] + perm[:a]


def inv_for_ab(perm, n, a, b):
    seq = [(v - b) % n for v in rotated(perm, a)]
    return inv_count(seq)


def full_min(perm, n):
    best = None
    for a in range(n):
        rot = rotated(perm, a)
        for b in range(n):
            seq = [(v - b) % n for v in rot]
            iv = inv_count(seq)
            if best is None or iv < best:
                best = iv
    return best


def fixed_a0_min(perm, n):
    best = None
    for b in range(n):
        seq = [(v - b) % n for v in perm]
        iv = inv_count(seq)
        if best is None or iv < best:
            best = iv
    return best


RULES = {
    "b=pi(a)": lambda perm, n, a: perm[a],
    "b=pi(a)-1": lambda perm, n, a: (perm[a] - 1) % n,
    "b=pi(a-1)+1": lambda perm, n, a: (perm[(a - 1) % n] + 1) % n,
}


def coupled_min(perm, n, rule):
    best = None
    for a in range(n):
        b = rule(perm, n, a)
        iv = inv_for_ab(perm, n, a, b)
        if best is None or iv < best:
            best = iv
    return best


def run(n, perms, exhaustive):
    bound = (n - 1) ** 2 // 4
    worst_full = 0
    worst_full_example = None
    worst_a0 = 0
    worst_a0_example = None
    worst_rule = {name: 0 for name in RULES}
    worst_rule_example = {name: None for name in RULES}
    checked = 0
    for perm in perms:
        p = list(perm)
        checked += 1
        if exhaustive:
            mf = full_min(p, n)
            if mf > worst_full:
                worst_full = mf
                worst_full_example = tuple(p)
        ma0 = fixed_a0_min(p, n)
        if ma0 > worst_a0:
            worst_a0 = ma0
            worst_a0_example = tuple(p)
        for name, rule in RULES.items():
            mr = coupled_min(p, n, rule)
            if mr > worst_rule[name]:
                worst_rule[name] = mr
                worst_rule_example[name] = tuple(p)
    return {
        "n": n,
        "bound": bound,
        "checked": checked,
        "exhaustive": exhaustive,
        "worst_full": worst_full if exhaustive else None,
        "worst_full_example": worst_full_example,
        "worst_fixed_a0": worst_a0,
        "worst_fixed_a0_example": worst_a0_example,
        "worst_rule": worst_rule,
        "worst_rule_example": worst_rule_example,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=None,
                     help="exhaustive over all permutations for 4..nmax")
    ap.add_argument("--n", type=int, default=None,
                     help="single n, exhaustive unless --sample given")
    ap.add_argument("--sample", type=int, default=None,
                     help="random sample size instead of exhaustive")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    results = []
    t0 = time.time()
    if args.nmax:
        for n in range(4, args.nmax + 1):
            perms = itertools.permutations(range(n))
            results.append(run(n, perms, exhaustive=True))
    elif args.n:
        if args.sample:
            rng = random.Random(args.seed)
            perms = []
            for _ in range(args.sample):
                p = list(range(args.n))
                rng.shuffle(p)
                perms.append(tuple(p))
            results.append(run(args.n, perms, exhaustive=False))
        else:
            perms = itertools.permutations(range(args.n))
            results.append(run(args.n, perms, exhaustive=True))
    else:
        raise SystemExit("need --nmax or --n")
    elapsed = time.time() - t0

    out_dir = os.path.join(ROOT, "data", "runs", "h13i_coupling")
    os.makedirs(out_dir, exist_ok=True)
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    lines = [f"\n## Run {ts} ({elapsed:.1f} s)\n"]
    lines.append(f"args: {vars(args)}\n")
    for r in results:
        cov = "exhaustive" if r["exhaustive"] else f"SAMPLED n={r['checked']}"
        lines.append(f"\n### n = {r['n']} ({cov}), bound floor((n-1)^2/4) = {r['bound']}\n")
        if r["exhaustive"]:
            lines.append(f"- worst I(pi) (full n^2 search) = {r['worst_full']}, "
                          f"example {r['worst_full_example']}\n")
        lines.append(f"- worst (fix a=0, vary b only) = {r['worst_fixed_a0']}, "
                      f"example {r['worst_fixed_a0_example']}\n")
        for name in RULES:
            lines.append(f"- worst (rule {name}) = {r['worst_rule'][name]}, "
                          f"example {r['worst_rule_example'][name]}\n")

    report_path = os.path.join(out_dir, "report.md")
    with open(report_path, "a") as f:
        f.writelines(lines)
    json_path = os.path.join(out_dir, f"result_{results[0]['n']}_{int(time.time())}.json")
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2)
    print("".join(lines))
    print(f"appended to {report_path}")


if __name__ == "__main__":
    main()
