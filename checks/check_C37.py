"""Checker for C37 (session 9): I(pi_h) <= floor((n-1)^2/4) on all reflections,
via an explicit two-block cut, for all n (not only 4 <= n <= 10).

  Part 1 (Lemma 0): block identity floor((n-1)^2/4) = C(floor(n/2),2) +
    C(ceil(n/2),2), exact integer arithmetic, 1 <= n <= 2000.
  Part 2 (Lemma 2): for pi_0(i) = (n-i) mod n, cut q = 0, shift s: direct
    inversion count of the relabelled line equals C(n-s,2) + C(s,2), for all
    s, 2 <= n <= 300 (brute-force inversion count, not the formula).
  Part 3 (Lemma 1): I(pi + t) = I(pi) for all t, checked by brute force over
    all n^2 cuts, random pi, 2 <= n <= 9.
  Part 4 (final bound): I(pi_h) <= floor((n-1)^2/4) for all n, all h, by full
    brute-force search over all n^2 cuts (not using lemmas 1-2 as a shortcut),
    2 <= n <= 10; where data/runs/line_profile/I_n{n}.bin is present (C33,
    4 <= n <= 10), also checks I(pi_h) == floor((n-1)^2/4) exactly (equality,
    not just <=).

Usage: python3 checks/check_C37.py [--nmax 10] [--random-cases 30]
Output: data/runs/check_C37/report.json, report.md.
"""

import argparse
import json
import os
import random
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
VERSION = "check_C37-1.0"
OUT = os.path.join(ROOT, "data", "runs", "check_C37")


def inv_count(w):
    n = len(w)
    c = 0
    for i in range(n):
        for j in range(i + 1, n):
            if w[i] > w[j]:
                c += 1
    return c


def line_at(pi, n, q, c):
    shift = (q + 1 - c) % n
    return [(pi[(q + 1 + j) % n] - shift) % n for j in range(n)]


def I_of(pi, n):
    best = None
    for q in range(n):
        for c in range(n):
            iv = inv_count(line_at(pi, n, q, c))
            if best is None or iv < best:
                best = iv
    return best


def target(n):
    return ((n - 1) ** 2) // 4


def block_identity(n):
    a = n // 2
    b = n - a
    lhs = target(n)
    rhs = a * (a - 1) // 2 + b * (b - 1) // 2
    return lhs, rhs


def lemma2_formula(n, s):
    return (n - s) * (n - s - 1) // 2 + s * (s - 1) // 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=10)
    ap.add_argument("--nmax-lemma0", type=int, default=2000)
    ap.add_argument("--nmax-lemma2", type=int, default=300)
    ap.add_argument("--random-cases", type=int, default=30)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    report = {"version": VERSION, "parts": {}}

    # Part 1: block identity
    fails0 = []
    for n in range(1, args.nmax_lemma0 + 1):
        lhs, rhs = block_identity(n)
        if lhs != rhs:
            fails0.append(n)
    report["parts"]["lemma0_block_identity"] = {
        "range": f"1<=n<={args.nmax_lemma0}",
        "fails": fails0,
        "pass": len(fails0) == 0,
    }
    print(f"Lemma 0: n=1..{args.nmax_lemma0}, fails={len(fails0)}")

    # Part 2: lemma 2 formula for pi_0, q=0
    fails2 = []
    for n in range(2, args.nmax_lemma2 + 1):
        pi0 = [(n - i) % n for i in range(n)]
        for s in range(n):
            w = line_at(pi0, n, 0, (1 - s) % n)
            iv = inv_count(w)
            expect = lemma2_formula(n, s)
            if iv != expect:
                fails2.append((n, s, iv, expect))
    report["parts"]["lemma2_formula"] = {
        "range": f"2<=n<={args.nmax_lemma2}, all s",
        "fails": fails2[:20],
        "n_fails": len(fails2),
        "pass": len(fails2) == 0,
    }
    print(f"Lemma 2: n=2..{args.nmax_lemma2}, all s, fails={len(fails2)}")

    # Part 3: I(pi+t) = I(pi), random pi, all t, brute force
    random.seed(args.seed)
    fails1 = []
    checked1 = 0
    for n in range(2, 10):
        for _ in range(args.random_cases):
            pi = list(range(n))
            random.shuffle(pi)
            Ibase = I_of(pi, n)
            for t in range(n):
                pit = [(v + t) % n for v in pi]
                It = I_of(pit, n)
                checked1 += 1
                if It != Ibase:
                    fails1.append((n, pi, t, Ibase, It))
    report["parts"]["lemma1_invariance"] = {
        "range": "2<=n<=9, random pi, all t",
        "checked": checked1,
        "fails": fails1[:20],
        "n_fails": len(fails1),
        "pass": len(fails1) == 0,
    }
    print(f"Lemma 1: checked={checked1}, fails={len(fails1)}")

    # Part 4: final bound on all reflections, all n, all h; plus equality vs C33 tables
    fails4 = []
    equality_checked = 0
    equality_mismatches = []
    for n in range(2, args.nmax + 1):
        tgt = target(n)
        ipath = os.path.join(ROOT, "data", "runs", "line_profile", f"I_n{n}.bin")
        Ivals = None
        if os.path.exists(ipath):
            Ivals = open(ipath, "rb").read()
        for h in range(n):
            pi = [(h - i) % n for i in range(n)]
            I = I_of(pi, n)
            if I > tgt:
                fails4.append((n, h, I, tgt))
            if Ivals is not None:
                # rank of pi in lexicographic order among permutations of range(n)
                import sys as _sys
                _sys.path.insert(0, os.path.join(ROOT, "exact"))
                from bfs import factorials, rank_perm  # noqa: E402
                fact = factorials(n)
                r = rank_perm(tuple(pi), fact)
                equality_checked += 1
                if Ivals[r] != I or Ivals[r] != tgt:
                    equality_mismatches.append((n, h, I, tgt, Ivals[r]))
    report["parts"]["final_bound"] = {
        "range": f"2<=n<={args.nmax}, all h",
        "fails": fails4,
        "pass": len(fails4) == 0,
        "equality_checked_against_line_profile": equality_checked,
        "equality_mismatches": equality_mismatches[:20],
        "equality_pass": len(equality_mismatches) == 0,
    }
    print(f"Final bound: n=2..{args.nmax}, all h, fails={len(fails4)}; "
          f"equality vs line_profile checked={equality_checked}, mismatches={len(equality_mismatches)}")

    report["elapsed_s"] = time.time() - t0
    overall = (
        report["parts"]["lemma0_block_identity"]["pass"]
        and report["parts"]["lemma2_formula"]["pass"]
        and report["parts"]["lemma1_invariance"]["pass"]
        and report["parts"]["final_bound"]["pass"]
        and report["parts"]["final_bound"]["equality_pass"]
    )
    report["overall_pass"] = overall

    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(report, f, indent=2, default=str)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write(f"# check_C37 report ({VERSION})\n\n")
        f.write(f"Overall: {'PASS' if overall else 'FAIL'}\n\n")
        for name, part in report["parts"].items():
            f.write(f"## {name}\n\n")
            f.write(f"range: {part.get('range')}\n\n")
            f.write(f"pass: {part.get('pass')}\n\n")
        f.write(f"\nElapsed: {report['elapsed_s']:.1f} s\n")

    print("OVERALL:", "PASS" if overall else "FAIL")
    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
