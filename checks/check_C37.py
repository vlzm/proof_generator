"""Checker for C37/C38 (session 9, H13-I attempt).

  C37 (arc/quadrant reformulation of double-cut inversions, PROVED): for a
    pair of positions i1 < i2 with delta = i2 - i1, v1 = pi(i1), v2 = pi(i2),
    eps = (v2 - v1) mod n, the pair is an inversion of the double rotation
    pi_{a,b}(j) = (pi((a + j) % n) - b) % n  (j = 0..n-1)
  iff exactly one of
    (a - i1) % n in {1, ..., delta}          ("i2 before i1" arc A2)
    (b - v1) % n in {1, ..., eps}             ("v1 > v2" arc B1)
  holds.  Part 1 checks this pair-by-pair against the brute-force inversion
  count of pi_{a,b}, for every (a, b), on random pi at each n.

  C38 (single-shift reductions of H13-I fail, VERIFIED counterexamples):
    fixing the position-rotation a by a cheap canonical rule (a = 0 always;
    a = position of value 0; a = that + 1; a = position of the median value
    (n-1)//2) and searching only the value-rotation b does not reach the
    bound floor((n-1)^2/4) for every pi -- explicit exhaustive counterexamples
    exist already at n = 7 (rule a = 0) and n = 8 (the other three rules).
    Part 2 re-derives these exhaustively and reports the counterexamples.

Usage: python3 checks/check_C37.py [--pairs_amax 9] [--rules_amax 8]
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
sys.path.insert(0, os.path.join(ROOT, "oracle"))
from moves import CORE_VERSION  # noqa: E402

VERSION = "check_C37-1.0"
OUT = os.path.join(ROOT, "data", "runs", "check_C37")


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def pi_ab(pi, a, b):
    n = len(pi)
    return tuple((pi[(a + j) % n] - b) % n for j in range(n))


def part1_pair_formula(amax, log, samples_per_n=20, seed=1):
    """Check the C37 arc/quadrant formula against brute force, all (a, b),
    for `samples_per_n` random permutations at each n (all n! for n <= 6)."""
    rng = random.Random(seed)
    ok = True
    rows = []
    for n in range(2, amax + 1):
        t0 = time.time()
        if n <= 6:
            pis = list(itertools.permutations(range(n)))
        else:
            pis = [tuple(rng.sample(range(n), n)) for _ in range(samples_per_n)]
        checked = 0
        for pi in pis:
            for a in range(n):
                for b in range(n):
                    w = pi_ab(pi, a, b)
                    inv_direct = inversions(w)
                    inv_formula = 0
                    for i1 in range(n):
                        for i2 in range(i1 + 1, n):
                            delta = i2 - i1
                            v1, v2 = pi[i1], pi[i2]
                            eps = (v2 - v1) % n
                            in_a2 = 1 <= (a - i1) % n <= delta
                            in_b1 = 1 <= (b - v1) % n <= eps
                            if in_a2 != in_b1:
                                inv_formula += 1
                    if inv_formula != inv_direct:
                        ok = False
                        log(f"FAIL C37 n={n} pi={pi} a={a} b={b}: formula={inv_formula} direct={inv_direct}")
            checked += 1
        rows.append({"n": n, "permutations_checked": checked, "all_a_b": True,
                     "coverage": "exhaustive" if n <= 6 else f"random samples_per_n={samples_per_n} seed={seed}",
                     "seconds": round(time.time() - t0, 2)})
        log(f"C37 n={n}: {checked} pi x {n * n} (a,b) pairs x C({n},2) inner pairs -- "
            f"{'exhaustive' if n <= 6 else 'sampled'}, {time.time() - t0:.2f} s")
    return ok, rows


def i_value_only(pi, a):
    n = len(pi)
    return min(inversions(pi_ab(pi, a, b)) for b in range(n))


RULES = {
    "a=0": lambda pi: 0,
    "a=pos(0)": lambda pi: pi.index(0),
    "a=pos(0)+1": lambda pi: (pi.index(0) + 1) % len(pi),
    "a=pos(median)": lambda pi: pi.index((len(pi) - 1) // 2),
}


def part2_single_shift_rules(amax, log):
    """For each canonical rule fixing a, find the exhaustive worst pi and
    check whether i_value_only(pi, rule(pi)) stays <= floor((n-1)^2/4)."""
    ok_overall = True
    rows = []
    for n in range(4, amax + 1):
        bound = (n - 1) ** 2 // 4
        t0 = time.time()
        worst = {name: (0, None) for name in RULES}
        for pi in itertools.permutations(range(n)):
            for name, rule in RULES.items():
                a = rule(pi)
                v = i_value_only(pi, a)
                if v > worst[name][0]:
                    worst[name] = (v, pi)
        row = {"n": n, "bound_floor_(n-1)^2/4": bound, "seconds": round(time.time() - t0, 1), "rules": {}}
        for name, (v, arg) in worst.items():
            fails = v > bound
            row["rules"][name] = {"max": v, "argmax": list(arg) if arg else None, "exceeds_bound": fails}
            log(f"C38 n={n} rule {name}: worst I_value_only = {v} vs bound {bound} at {arg} "
                f"-> {'EXCEEDS (rule insufficient)' if fails else 'within bound'}")
        rows.append(row)
    # C38 claim: each rule has a counterexample by n = amax (already found at n <= 8)
    any_exceeded = {name: any(r["rules"][name]["exceeds_bound"] for r in rows) for name in RULES}
    log(f"C38 summary: rules with a confirmed counterexample (n <= {amax}): "
        f"{[k for k, v in any_exceeded.items() if v]}")
    return rows, any_exceeded


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs_amax", type=int, default=9)
    ap.add_argument("--rules_amax", type=int, default=8)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    lines = []

    def log(msg):
        print(msg, flush=True)
        lines.append(msg)

    log(f"== {VERSION} {CORE_VERSION} args={vars(args)}")
    ok1, rows1 = part1_pair_formula(args.pairs_amax, log)
    rows2, any_exceeded = part2_single_shift_rules(args.rules_amax, log)
    verdict = "PASS" if ok1 else "FAIL"
    log(f"verdict (C37 formula correctness): {verdict}")
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump({"version": VERSION, "core": CORE_VERSION, "args": vars(args),
                    "C37_pair_formula": rows1, "C37_ok": ok1,
                    "C38_single_shift_rules": rows2, "C38_any_counterexample": any_exceeded,
                    "verdict": verdict}, f, indent=1)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# check_C37/C38 -- arc/quadrant reformulation, single-shift rules ruled out\n\n")
        f.write(f"Команда: `python3 checks/check_C37.py --pairs_amax {args.pairs_amax} "
                f"--rules_amax {args.rules_amax}`. Версии: {VERSION}, {CORE_VERSION}.\n\n```text\n")
        f.write("\n".join(lines) + "\n```\n")
    sys.exit(0 if verdict == "PASS" else 1)


if __name__ == "__main__":
    main()
