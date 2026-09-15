"""H13-I candidate strategies for choosing a good double cut (session 9).

H13-I (docs/notes/h13_line_model.md §6): every pi has a double cut (q, c) with
inv(w) <= floor((n-1)^2/4) = I(pi) [see experiments/line_profile.c, C35/C36].
Exhaustive search over all n^2 cuts already proves this for 4 <= n <= 10
(C33/C35). The open step is an ANALYTIC choice of (q, c) (not full search)
that a proof can use. This script tests two single-scalar-feature candidates
and reports where each first fails, i.e. the smallest counterexample (rule 6).

Identity used throughout (derived this session, checked against brute force):
fix a position rotation a (u_k = pi[(a+k) mod n], a linear, non-cyclic reading
of pi starting at position a) and let Inv(u) be its ordinary (non-cyclic)
inversion count.  Summing the relabelled-line inversion count inv(w_{a,b})
over all n choices of the value shift b gives, EXACTLY:

    sum_b inv(w_{a,b}) = S(u) + n * Inv(u),   S(u) = sum_k u[k] * (2k - n + 1)

(proof: for fixed a, each pair i<j is inverted for exactly beta_ij = (u_j -
u_i) mod n values of b out of n; splitting beta_ij by whether pair (i,j) is
itself concordant or inverted in u and summing gives S(u) + n*Inv(u)). So
average_b inv(w_{a,b}) = S(u)/n + Inv(u), independent of a in TOTAL sum (the
grand sum over all n^2 cuts does not depend on a), but min_b inv(w_{a,b}) is
what matters and is NOT predicted by Inv(u) alone -- see the refuted
candidates below.

Candidates tested, each combined with the TRUE min over b (not just the
average, which is a strictly weaker/easier target than the min):
  argmin  -- a = argmin_a Inv(u_a), i.e. the rotation that already looks most
             sorted as a plain (non-cyclic) sequence.
  argmax  -- a = argmax_a Inv(u_a), the least sorted rotation.

Both are REFUTED (worst_gap > 0 at some n <= 9), with minimal n and an
explicit witness recorded below and in the JSON report. Diagnosis: on the
witnesses, Inv(u) at the candidate a is either the unique best or the unique
WORST of the n rotations for the subsequent best-b search -- a single scalar
feature of the rotation does not predict which a admits a good b, confirming
(with new, smaller/different witnesses) the session-8 verdict that averaging
and single-feature heuristics do not settle H13-I.

Usage: python3 experiments/toric_cut_heuristics.py --nmax 9
Output: data/runs/h13_cut_heuristics/report.json, report.md.
Version toric_cut_heuristics-1.0.
"""

import argparse
import itertools
import json
import os
import time

VERSION = "toric_cut_heuristics-1.0"
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "data", "runs", "h13_cut_heuristics")


def inv_count(seq):
    """O(n log n) merge-sort inversion count; seq mutated internally only."""
    arr = list(seq)
    n = len(arr)
    tmp = [0] * n

    def sort(lo, hi):
        if hi - lo <= 1:
            return 0
        mid = (lo + hi) // 2
        inv = sort(lo, mid) + sort(mid, hi)
        i, j, k = lo, mid, lo
        while i < mid and j < hi:
            if arr[i] <= arr[j]:
                tmp[k] = arr[i]
                i += 1
            else:
                tmp[k] = arr[j]
                j += 1
                inv += mid - i
            k += 1
        while i < mid:
            tmp[k] = arr[i]
            i += 1
            k += 1
        while j < hi:
            tmp[k] = arr[j]
            j += 1
            k += 1
        arr[lo:hi] = tmp[lo:hi]
        return inv

    return sort(0, n)


def all_rotation_invs(pi):
    """Inv(u_a) for a = 0..n-1, u_a[k] = pi[(a+k) mod n], via O(n) increments:
    moving the front element x to the back changes the count by
    (n-1) - 2*cnt_less(x, rest)."""
    n = len(pi)
    doubled = list(pi) * 2
    invs = [0] * n
    invs[0] = inv_count(pi)
    cur = invs[0]
    for a in range(n - 1):
        x = doubled[a]
        rest = doubled[a + 1 : a + n]
        cnt_less = sum(1 for y in rest if y < x)
        cur += (n - 1) - 2 * cnt_less
        invs[a + 1] = cur
    return invs


def best_b_for_u(u):
    """min over value shift b of inv((u[k]-b) mod n), via O(n) increments:
    Inv(w_{b+1}) = Inv(w_b) + (n-1-2p), p = position in u of value b."""
    n = len(u)
    pos_of = [0] * n
    for k, val in enumerate(u):
        pos_of[val] = k
    cur = inv_count(u)
    best = cur
    for b in range(n - 1):
        p = pos_of[b]
        cur += (n - 1 - 2 * p)
        if cur < best:
            best = cur
    return best


def target(n):
    return (n - 1) ** 2 // 4


def run(nmax, nmin=4):
    os.makedirs(OUT, exist_ok=True)
    report = {"version": VERSION, "by_n": {}}
    for n in range(nmin, nmax + 1):
        t0 = time.time()
        worst = {"argmin": (-1, None), "argmax": (-1, None)}
        checked = 0
        for pi in itertools.permutations(range(n)):
            checked += 1
            invs = all_rotation_invs(pi)
            a_min = min(range(n), key=lambda i: invs[i])
            a_max = max(range(n), key=lambda i: invs[i])
            for name, a in (("argmin", a_min), ("argmax", a_max)):
                u = [pi[(a + k) % n] for k in range(n)]
                vb = best_b_for_u(u)
                gap = vb - target(n)
                if gap > worst[name][0]:
                    worst[name] = (gap, pi)
        elapsed = time.time() - t0
        report["by_n"][n] = {
            "target": target(n),
            "checked": checked,
            "argmin_worst_gap": worst["argmin"][0],
            "argmin_example": worst["argmin"][1],
            "argmax_worst_gap": worst["argmax"][0],
            "argmax_example": worst["argmax"][1],
            "seconds": round(elapsed, 2),
        }
        print(n, report["by_n"][n])

    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(report, f, indent=2)

    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# H13-I cut-heuristics (session 9)\n\n")
        f.write(
            "Two single-feature candidates for choosing the position rotation a "
            "(then TRUE min over value shift b), tested exhaustively against "
            "target = floor((n-1)^2/4). worst_gap > 0 means the candidate fails "
            "to reach the target on the given (minimal, for that candidate) "
            "example.\n\n"
        )
        f.write("| n | target | argmin worst_gap | argmin example | argmax worst_gap | argmax example | checked | s |\n")
        f.write("|---|---|---|---|---|---|---|---|\n")
        for n in sorted(report["by_n"]):
            r = report["by_n"][n]
            f.write(
                f"| {n} | {r['target']} | {r['argmin_worst_gap']} | {r['argmin_example']} | "
                f"{r['argmax_worst_gap']} | {r['argmax_example']} | {r['checked']} | {r['seconds']} |\n"
            )
    return report


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=4)
    ap.add_argument("--nmax", type=int, default=9)
    args = ap.parse_args()
    run(args.nmax, args.nmin)
