"""Sampled screening of H13-I for large n (session 9, PLAN.md §0/§8).

H13-I: for every permutation pi of Z_n there is a double cut (q, c) whose
relabelled line has at most floor((n-1)^2/4) inversions, i.e.
I(pi) = min_{q,c} inv(w) <= floor((n-1)^2/4).  Exhaustive verification exists
only for 4 <= n <= 10 (experiments/line_profile.c, C33/C35).  Exhaustive
enumeration is infeasible for n = 11 (see PLAN.md, "n = 11 не запускается
как следующая строка"); this script instead follows AGENTS.md rule 9: sampled
and structured inputs at n = 20, 50, 100, 101 and around thresholds, up to
n = 200, looking for a violation (a bug hunt, not a proof).

I(pi) is computed in O(n^2) per permutation (O(n) rotations of the position
cut, O(n) incremental update of inversions over the value cut within each
rotation, after one O(n log n) merge-sort inversion count) instead of the
O(n^4) brute force in line_profile.c; verified to agree with brute force
(hence with line_profile.c) for all permutations, 2 <= n <= 8.

Usage: python3 experiments/i_bound_sampled.py [--nmax 200] [--randoms 30] [--seed 12345]
Output: data/runs/i_bound_sampled/report.json, report.md.  Version 1.0.
"""

import argparse
import json
import os
import random
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
VERSION = "i_bound_sampled-1.0"
OUT = os.path.join(ROOT, "data", "runs", "i_bound_sampled")


def merge_inversions(v):
    n = len(v)
    arr = list(v)
    tmp = [0] * n
    inv = [0]

    def rec(lo, hi):
        if hi - lo <= 1:
            return
        mid = (lo + hi) // 2
        rec(lo, mid)
        rec(mid, hi)
        i, j, k = lo, mid, lo
        cnt = 0
        while i < mid and j < hi:
            if arr[i] <= arr[j]:
                tmp[k] = arr[i]; i += 1
            else:
                tmp[k] = arr[j]; j += 1
                cnt += mid - i
            k += 1
        while i < mid:
            tmp[k] = arr[i]; i += 1; k += 1
        while j < hi:
            tmp[k] = arr[j]; j += 1; k += 1
        arr[lo:hi] = tmp[lo:hi]
        inv[0] += cnt

    rec(0, n)
    return inv[0]


def I_of_pi(pi):
    """min over (a, b) in Z_n x Z_n of inv(w), w_j = pi[(a+j) % n] - b mod n."""
    n = len(pi)
    best = None
    for a in range(n):
        v = [pi[(a + j) % n] for j in range(n)]
        pos = [0] * n
        for idx, val in enumerate(v):
            pos[val] = idx
        cur = merge_inversions(v)
        if best is None or cur < best:
            best = cur
        for b in range(1, n):
            p = pos[b - 1]
            cur += n - 1 - 2 * p
            if cur < best:
                best = cur
    return best


def target(n):
    return (n - 1) * (n - 1) // 4


def build_inputs(n, randoms, seed):
    rnd = random.Random(seed * 100000 + n)
    inputs = []
    sigma = [(1 - i) % n for i in range(n)]
    inputs.append(("sigma_n", sigma))
    for k in (1, 2, n // 3, n // 2):
        inputs.append((f"rot_sigma_k{k}", [sigma[(i + k) % n] for i in range(n)]))
    a_candidates = sorted(set([2, 3, 5, n - 2, n - 3, (n // 2) + (1 if n % 2 == 0 else 0)]))
    for a in a_candidates:
        for b in (0, 1, n // 4, n // 2, n - 1):
            pi = [(a * i + b) % n for i in range(n)]
            if len(set(pi)) == n:
                inputs.append((f"affine_a{a}_b{b}", pi))
    for t in range(randoms):
        pi = list(range(n))
        rnd.shuffle(pi)
        inputs.append((f"random{t}", pi))
    inputs.append(("rev_n", list(reversed(range(n)))))
    return inputs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nlist", type=str,
                     default="11,12,13,14,15,16,17,18,19,20,21,22,24,27,30,40,50,64,80,100,101,127,150,200")
    ap.add_argument("--randoms", type=int, default=30)
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()
    ns = [int(x) for x in args.nlist.split(",") if x]

    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    total_checks = 0
    violations = []
    by_n = []
    for n in ns:
        best_gap = None
        tight_labels = []
        for label, pi in build_inputs(n, args.randoms, args.seed):
            I = I_of_pi(pi)
            t = target(n)
            total_checks += 1
            gap = I - t
            if best_gap is None or gap > best_gap:
                best_gap = gap
                tight_labels = [label]
            elif gap == best_gap:
                tight_labels.append(label)
            if I > t:
                violations.append({"n": n, "label": label, "I": I, "target": t})
        by_n.append({"n": n, "max_I_minus_target": best_gap, "tight_at": tight_labels[:5]})

    report = {
        "version": VERSION,
        "ns": ns,
        "randoms_per_n": args.randoms,
        "seed": args.seed,
        "total_checks": total_checks,
        "violations": violations,
        "by_n": by_n,
        "elapsed_s": time.time() - t0,
    }
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(report, f, indent=2)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write(f"# {VERSION} report\n\n")
        f.write(f"H13-I sampled screening (rule 9): random, affine, sigma_n + rotations, "
                f"reversal at n in {ns}, {args.randoms} random inputs per n.\n\n")
        f.write(f"Total checks: {total_checks}. Violations: {len(violations)}.\n\n")
        f.write("| n | max(I - target) among sampled | tight at |\n|---|---|---|\n")
        for row in by_n:
            f.write(f"| {row['n']} | {row['max_I_minus_target']} | {', '.join(row['tight_at'])} |\n")
        if violations:
            f.write("\n## Violations\n\n")
            for v in violations:
                f.write(f"- n={v['n']} {v['label']}: I={v['I']} > target={v['target']}\n")
    print(f"total_checks={total_checks} violations={len(violations)} elapsed={time.time()-t0:.1f}s")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
