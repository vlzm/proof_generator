"""Checker for C37 (pair decomposition of the toric inversion sum, session 9).

Claim (docs/proofs/C37_pair_decomposition.md): for a permutation pi of
{0,...,n-1}, define for the double cut (q, c) the line
w_j = pi[(q+1+j) % n] - (q+1-c) mod n (PROBLEM Sec. 5.6). Then

    S(pi) := sum over all n^2 cuts (q, c) of inv(w)
           = sum_{0<=i<j<=n-1} F(j-i, (pi[j]-pi[i]) mod n),
    F(d, e) = n*(d+e) - 2*d*e.

This is checked two ways:
  (a) exhaustively over all permutations for 4 <= n <= AMAX, comparing the
      direct O(n^4) computation of S(pi) against the closed-form pair sum;
  (b) on random permutations for larger n (up to RMAX), same comparison.

The checker also re-derives, per pair, the number of q with "i before j" and
the number of c with "pi(i) before pi(j) in rank order", and checks these
equal n - d and n - e exactly (the two counting lemmas used in the proof),
on all pairs of a few sample permutations.

This does NOT prove or refute H13-I; see docs/notes/h13_line_model.md Sec. 7
for the verdict of session 9.

Usage: python3 checks/check_C37.py [--amax 8] [--rmax 60] [--samples 20]
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


def S_direct(pi):
    n = len(pi)
    total = 0
    for q in range(n):
        for c in range(n):
            shift = (q + 1 - c) % n
            w = [(pi[(q + 1 + j) % n] - shift) % n for j in range(n)]
            total += inversions(w)
    return total


def F(d, e, n):
    return n * (d + e) - 2 * d * e


def S_formula(pi):
    n = len(pi)
    total = 0
    for i in range(n):
        for j in range(i + 1, n):
            d = j - i
            e = (pi[j] - pi[i]) % n
            total += F(d, e, n)
    return total


def count_bits(pi):
    """Per-pair check of the two counting lemmas: #{q: P=1} = n-d, #{c: Q=1} = n-e."""
    n = len(pi)
    ok = True
    for i in range(n):
        for j in range(i + 1, n):
            d = j - i
            e = (pi[j] - pi[i]) % n
            # P(t): i before j at line-index (p - t) mod n, t = q+1
            cntP = sum(1 for t in range(n) if (i - t) % n < (j - t) % n)
            # Q(c): rank_c(pi[i]) < rank_c(pi[j])
            cntQ = sum(1 for c in range(n) if (pi[i] - c) % n < (pi[j] - c) % n)
            if cntP != n - d or cntQ != n - e:
                ok = False
    return ok


def part_a(amax, log):
    rows = []
    ok = True
    for n in range(4, amax + 1):
        t0 = time.time()
        cnt = 0
        for pi in itertools.permutations(range(n)):
            s1 = S_direct(pi)
            s2 = S_formula(pi)
            if s1 != s2:
                ok = False
                log(f"FAIL part A: n={n} pi={pi} S_direct={s1} S_formula={s2}")
                break
            cnt += 1
        rows.append({"n": n, "count": cnt, "seconds": round(time.time() - t0, 1)})
        log(f"part A n={n}: {cnt} permutations, S_direct == S_formula, {time.time() - t0:.1f} s")
    return ok, rows


def part_b(rmax, samples, log, seed=1):
    rng = random.Random(seed)
    rows = []
    ok = True
    n = 9
    while n <= rmax:
        for _ in range(samples):
            pi = list(range(n))
            rng.shuffle(pi)
            pi = tuple(pi)
            s1 = S_direct(pi)
            s2 = S_formula(pi)
            if s1 != s2:
                ok = False
                log(f"FAIL part B: n={n} pi={pi} S_direct={s1} S_formula={s2}")
        rows.append({"n": n, "samples": samples})
        log(f"part B n={n}: {samples} random permutations, S_direct == S_formula")
        n += 1 if n < 20 else (10 if n < 40 else 20)
    return ok, rows


def part_c(log):
    """Per-pair counting lemma (#P=1 is n-d, #Q=1 is n-e), exhaustively for small n."""
    ok = True
    checked = 0
    for n in range(4, 7):
        for pi in itertools.permutations(range(n)):
            if not count_bits(pi):
                ok = False
                log(f"FAIL part C: n={n} pi={pi}")
            checked += 1
    log(f"part C: per-pair counting lemma checked on all permutations 4 <= n <= 6 ({checked} total)")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amax", type=int, default=8)
    ap.add_argument("--rmax", type=int, default=60)
    ap.add_argument("--samples", type=int, default=20)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    lines = []

    def log(msg):
        print(msg, flush=True)
        lines.append(msg)

    log(f"== {VERSION} {CORE_VERSION} args={vars(args)}")
    oka, ra = part_a(args.amax, log)
    okb, rb = part_b(args.rmax, args.samples, log)
    okc = part_c(log)
    verdict = "PASS" if (oka and okb and okc) else "FAIL"
    log(f"verdict: {verdict}")

    report = {"version": VERSION, "core": CORE_VERSION, "args": vars(args),
              "part_a": ra, "part_b": rb, "part_c_ok": okc, "verdict": verdict}
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(report, f, indent=2)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# check_C37 report\n\n```\n" + "\n".join(lines) + "\n```\n")
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
