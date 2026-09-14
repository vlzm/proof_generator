"""Checker for C37 (session 9, H13-I exploration): the cut-exchange recurrence
  for the double-cut concordant-pair count, and the insufficiency of two
  restricted one-parameter search families.

  Setup (docs/notes/h13_line_model.md Section 0): pi a permutation of Z_n,
  X = {(i, pi(i)) : i in Z_n} its graph on the torus Z_n x Z_n.  For a
  position cut a and a value cut b, linearize positions starting at a
  (p(i) = (i - a) mod n) and values starting at b (v(i) = (pi(i) - b) mod n).
  A pair of points i != i' is *concordant* at (a, b) iff p(i) < p(i') and
  v(i) < v(i'), or p(i) > p(i') and v(i) > v(i').  F(a, b) = number of
  concordant pairs.  Writing w_j = pi((a + j) mod n) - b mod n for the
  resulting linear sequence, F(a, b) = C(n, 2) - inv(w).  This is exactly the
  (q, c)-double-cut of the line model with a = q + 1, b = q + 1 - c.

  Part A -- recurrence (elementary, proved in docs/notes/h13_line_model.md
    session 9 addendum by a direct exchange argument, no computation
    required for correctness; verified here against brute force):

      F(a+1, b) - F(a, b) = 2 * V(a, b) - (n - 1),   V(a, b) = (pi(a) - b) mod n
      F(a, b+1) - F(a, b) = 2 * U(a, b) - (n - 1),   U(a, b) = (pi^-1(b) - a) mod n

  Part B -- restricted search families that do NOT always reach the H13-I
    target floor((n-1)^2/4) (exhaustive over all pi, small n):
      (b1) a = 0 fixed, vary b only (n candidates);
      (b2) diagonal cuts b = pi(a) for each a (n candidates, cut placed at a
           data point).
    Both fail on n = 7 (b1 already fails on the Lemma C example from the
    note, (0,5,4,3,2,1,6)); this documents why the full n^2 search over
    (a, b) cannot be collapsed to either restricted family.

  Part C -- sanity re-check of the already-registered C33/C35 fact that the
    full 2D search does reach the target for 4 <= n <= AMAX (no new claim;
    cross-check only).

Usage: python3 checks/check_C37.py [--amax 7] [--random-nmax 12] [--trials 50]
Output: data/runs/check_C37/report.json, report.md.
"""

import argparse
import itertools
import json
import os
import random
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "data", "runs", "check_C37")

VERSION = "check_C37-1.0"


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def linearize(pi, a, b):
    n = len(pi)
    return [(pi[(a + j) % n] - b) % n for j in range(n)]


def F(pi, a, b):
    n = len(pi)
    total = n * (n - 1) // 2
    return total - inversions(linearize(pi, a, b))


def inv_of(pi):
    n = len(pi)
    inv = [0] * n
    for i, v in enumerate(pi):
        inv[v] = i
    return inv


def part_a_recurrence(amax, random_nmax, trials, log):
    """Verify F(a+1,b)-F(a,b) and F(a,b+1)-F(a,b) against brute force."""
    ok = True
    checked_pairs = 0
    tested_ns = []

    def check_one(pi):
        nonlocal ok, checked_pairs
        n = len(pi)
        piinv = inv_of(pi)
        for a in range(n):
            for b in range(n):
                lhs1 = F(pi, (a + 1) % n, b) - F(pi, a, b)
                V = (pi[a] - b) % n
                rhs1 = 2 * V - (n - 1)
                lhs2 = F(pi, a, (b + 1) % n) - F(pi, a, b)
                U = (piinv[b] - a) % n
                rhs2 = 2 * U - (n - 1)
                checked_pairs += 2
                if lhs1 != rhs1 or lhs2 != rhs2:
                    ok = False
                    log.append(f"MISMATCH pi={pi} a={a} b={b} "
                               f"lhs1={lhs1} rhs1={rhs1} lhs2={lhs2} rhs2={rhs2}")

    for n in range(4, amax + 1):
        tested_ns.append(("exhaustive", n))
        for pi in itertools.permutations(range(n)):
            check_one(list(pi))

    rng = random.Random(20260913)
    for n in range(amax + 1, random_nmax + 1):
        tested_ns.append(("random", n, trials))
        for _ in range(trials):
            pi = list(range(n))
            rng.shuffle(pi)
            check_one(pi)

    return {
        "ok": ok,
        "checked_pairs": checked_pairs,
        "coverage": tested_ns,
    }


def min_inv_a0(pi):
    n = len(pi)
    return min(inversions(linearize(pi, 0, b)) for b in range(n))


def min_inv_diag(pi):
    n = len(pi)
    return min(inversions(linearize(pi, a, pi[a])) for a in range(n))


def min_inv_full(pi):
    n = len(pi)
    return min(inversions(linearize(pi, a, b)) for a in range(n) for b in range(n))


def part_b_restricted_families(amax, log):
    rows = []
    for n in range(4, amax + 1):
        target = ((n - 1) ** 2) // 4
        worst_a0 = -1
        worst_a0_perm = None
        exceed_a0 = 0
        worst_diag = -1
        worst_diag_perm = None
        exceed_diag = 0
        count = 0
        for pi in itertools.permutations(range(n)):
            pi = list(pi)
            count += 1
            v0 = min_inv_a0(pi)
            if v0 > target:
                exceed_a0 += 1
            if v0 > worst_a0:
                worst_a0, worst_a0_perm = v0, tuple(pi)
            vd = min_inv_diag(pi)
            if vd > target:
                exceed_diag += 1
            if vd > worst_diag:
                worst_diag, worst_diag_perm = vd, tuple(pi)
        rows.append({
            "n": n,
            "target": target,
            "total_perms": count,
            "a0_only": {
                "worst": worst_a0,
                "worst_perm": worst_a0_perm,
                "exceeding": exceed_a0,
            },
            "diagonal_only": {
                "worst": worst_diag,
                "worst_perm": worst_diag_perm,
                "exceeding": exceed_diag,
            },
        })
        log.append(
            f"n={n} target={target}: a=0-only worst={worst_a0} "
            f"(exceeds on {exceed_a0}/{count}); diagonal-only worst={worst_diag} "
            f"(exceeds on {exceed_diag}/{count})"
        )
    return rows


def part_c_full_search_sanity(amax, log):
    rows = []
    for n in range(4, amax + 1):
        target = ((n - 1) ** 2) // 4
        worst = -1
        for pi in itertools.permutations(range(n)):
            v = min_inv_full(list(pi))
            if v > worst:
                worst = v
        ok = worst <= target
        rows.append({"n": n, "target": target, "worst_full_search": worst, "ok": ok})
        log.append(f"n={n}: full 2D search worst min-inv = {worst} (target {target}) -> {'OK' if ok else 'FAIL'}")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amax", type=int, default=7, help="exhaustive n upper bound for parts A/B/C")
    ap.add_argument("--random-nmax", type=int, default=12, help="upper n for random spot checks in part A")
    ap.add_argument("--trials", type=int, default=50, help="random trials per n in part A beyond amax")
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    log = []
    t0 = time.time()

    log.append(f"{VERSION} amax={args.amax} random_nmax={args.random_nmax} trials={args.trials}")

    resA = part_a_recurrence(args.amax, args.random_nmax, args.trials, log)
    resB = part_b_restricted_families(args.amax, log)
    resC = part_c_full_search_sanity(args.amax, log)

    elapsed = time.time() - t0
    all_ok = resA["ok"] and all(r["ok"] for r in resC)

    report = {
        "version": VERSION,
        "elapsed_sec": elapsed,
        "part_a_recurrence": resA,
        "part_b_restricted_families": resB,
        "part_c_full_search_sanity": resC,
        "overall": "PASS" if all_ok else "FAIL",
    }

    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(report, f, indent=2, default=str)

    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write(f"# {VERSION} report\n\n")
        f.write(f"Overall: **{report['overall']}**  (elapsed {elapsed:.1f} s)\n\n")
        f.write("## Part A -- cut-exchange recurrence\n\n")
        f.write(f"Checked {resA['checked_pairs']} (a,b) recurrence instances, "
                f"coverage: {resA['coverage']}. Result: {'PASS' if resA['ok'] else 'FAIL'}\n\n")
        f.write("## Part B -- restricted one-parameter families (a=0 only; diagonal b=pi(a))\n\n")
        f.write("| n | target | a=0-only worst | a=0-only #exceeding | diagonal worst | diagonal #exceeding |\n")
        f.write("|---|---|---|---|---|---|\n")
        for r in resB:
            f.write(f"| {r['n']} | {r['target']} | {r['a0_only']['worst']} | "
                    f"{r['a0_only']['exceeding']}/{r['total_perms']} | "
                    f"{r['diagonal_only']['worst']} | "
                    f"{r['diagonal_only']['exceeding']}/{r['total_perms']} |\n")
        f.write("\n## Part C -- full 2D search sanity re-check (not a new claim; cross-check of C33)\n\n")
        f.write("| n | target | worst over all pi of min_{a,b} inv | OK |\n")
        f.write("|---|---|---|---|\n")
        for r in resC:
            f.write(f"| {r['n']} | {r['target']} | {r['worst_full_search']} | {r['ok']} |\n")

    for line in log:
        print(line)
    print(f"\nOverall: {report['overall']}  ({elapsed:.1f} s)")
    print(f"Report written to {OUT}/report.md")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
