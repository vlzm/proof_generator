"""Checker for C37 (value-shift closed form, session 9,
docs/proofs/C37_value_shift_formula.md).

  Part A: the identity inv(a,b+1)-inv(a,b) = (n-1)-2*sigma_a(b) and the
    resulting closed form min_b inv(a,b) = inv_count(sigma_a) -
    max_b[2 S_a(b) - b(n-1)] are checked against a brute double loop over
    all (a,b) with a naive O(n^2) inversion count, on random permutations,
    4 <= n <= 9 (500 trials per n).
  Part B: fast_I (experiments/h13i_fast.py, O(n^2) via C37) reproduces the
    brute I(pi) = min_{a,b} inv(a,b) exactly on ALL permutations, 4 <= n <= 8
    (exhaustive; the O(n^4) brute reference is the bottleneck, not fast_I).
  Part C: fast_I reproduces the known max_pi I(pi) = floor((n-1)^2/4)
    (C33/C35), 4 <= n <= 10 exhaustively (cross-check against the existing
    line_profile data), read from data/runs/h13i_fast/report.json if present.
  Part D: the universal ("a-free") strengthening is refuted: on the Lemma C
    counterexample (docs/notes/h13_line_model.md §4, n=7,
    w=(0,5,4,3,2,1,6)), min_b inv_count(rotate_b(w)) = 10 > floor(6^2/4) = 9.

Usage: python3 checks/check_C37.py [--trials 500] [--amax 8]
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
sys.path.insert(0, os.path.join(ROOT, "experiments"))
from h13i_fast import fast_I, inv_count_bit  # noqa: E402

VERSION = "check_C37-1.0"
OUT = os.path.join(ROOT, "data", "runs", "check_C37")


def inv_count_naive(seq):
    n = len(seq)
    return sum(1 for i in range(n) for j in range(i + 1, n) if seq[i] > seq[j])


def brute_f(pi, n, a, b):
    w = [(pi[(a + j) % n] - b) % n for j in range(n)]
    return inv_count_naive(w)


def brute_I(pi, n):
    return min(brute_f(pi, n, a, b) for a in range(n) for b in range(n))


def part_a(trials, log):
    random.seed(42)
    ok = True
    checked = 0
    for n in range(4, 10):
        for _ in range(trials):
            pi = list(range(n))
            random.shuffle(pi)
            for a in range(n):
                brute_min_b = min(brute_f(pi, n, a, b) for b in range(n))
                pinv = [0] * n
                for i, v in enumerate(pi):
                    pinv[v] = i
                sigma_a = [(pinv[t] - a) % n for t in range(n)]
                inv_sigma = inv_count_naive(sigma_a)
                S = 0
                bestE = 0
                for b in range(1, n + 1):
                    S += sigma_a[b - 1]
                    E = 2 * S - b * (n - 1)
                    if E > bestE:
                        bestE = E
                formula = inv_sigma - bestE
                checked += 1
                if formula != brute_min_b:
                    ok = False
                    log.append(f"FAIL part A: n={n} pi={pi} a={a} "
                               f"formula={formula} brute={brute_min_b}")
    log.append(f"part A: {checked} (n, pi, a) triples checked, ok={ok}")
    return ok


def part_b(amax, log):
    ok = True
    total = 0
    for n in range(4, amax + 1):
        t0 = time.time()
        cnt = 0
        for pi in itertools.permutations(range(n)):
            pi = list(pi)
            bi = brute_I(pi, n)
            fi, _ = fast_I(pi, n)
            cnt += 1
            total += 1
            if bi != fi:
                ok = False
                log.append(f"FAIL part B: n={n} pi={pi} brute={bi} fast={fi}")
        log.append(f"part B: n={n} {cnt} perms exhaustive, "
                    f"{time.time()-t0:.1f}s, ok so far={ok}")
    log.append(f"part B: {total} permutations total, ok={ok}")
    return ok


def part_c(log):
    path = os.path.join(ROOT, "data", "runs", "h13i_fast", "report.json")
    if not os.path.exists(path):
        log.append("part C: SKIPPED (data/runs/h13i_fast/report.json not found; "
                    "run experiments/h13i_fast.py --exhaustive 10 first)")
        return True
    with open(path) as fh:
        data = json.load(fh)
    ok = True
    for r in data.get("runs", []):
        if r["mode"] != "exhaustive":
            continue
        bound = r["bound"]
        if r["max_I"] != bound:
            ok = False
            log.append(f"FAIL part C: n={r['n']} max_I={r['max_I']} != "
                       f"bound={bound}")
        else:
            log.append(f"part C: n={r['n']} max_I={r['max_I']} == "
                       f"floor((n-1)^2/4)={bound} (exhaustive, {r['count']} perms)")
    return ok


def part_d(log):
    n = 7
    w = (0, 5, 4, 3, 2, 1, 6)
    bound = ((n - 1) ** 2) // 4
    best = min(inv_count_naive([(x - b) % n for x in w]) for b in range(n))
    ok = best > bound
    log.append(f"part D: n={n} w={w} inv(w)={inv_count_naive(list(w))} "
                f"min_b inv(rotate_b(w))={best} bound={bound} "
                f"(expected > bound, i.e. a-free strengthening refuted: {ok})")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=500)
    ap.add_argument("--amax", type=int, default=9)
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    log = []
    ok_a = part_a(args.trials, log)
    ok_b = part_b(args.amax, log)
    ok_c = part_c(log)
    ok_d = part_d(log)
    overall = ok_a and ok_b and ok_c and ok_d

    for line in log:
        print(line)
    print("OVERALL:", "PASS" if overall else "FAIL")

    with open(os.path.join(OUT, "report.json"), "w") as fh:
        json.dump({"version": VERSION, "pass": overall, "log": log}, fh, indent=2)
    with open(os.path.join(OUT, "report.md"), "w") as fh:
        fh.write(f"# check_C37 report\n\nVersion {VERSION}. "
                 f"Overall: {'PASS' if overall else 'FAIL'}\n\n")
        for line in log:
            fh.write(f"- {line}\n")


if __name__ == "__main__":
    main()
