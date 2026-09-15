"""Checker for C37/C38 (session 9, H13-I attempt): cutoff-shift calculus for
the double-cut line model of docs/notes/h13_line_model.md.

Definitions (PLAN H13, docs/notes/h13_line_model.md section 0): for pi in S_n,
edge cut q in Z_n and value shift c in Z_n, the line is
    w_j = (pi((q + 1 + j) mod n) - (q + 1 - c)) mod n,   j = 0 .. n-1.
inv(w) is the number of inversions of w read as a linear sequence.

Part A (C37): two exact identities for how inv(w) changes under a unit shift
of c or q (proved in docs/proofs/H13I_cutoff_calculus.md, Lemmas 1-2), plus the
inverse-duality identity (Lemma 3): inv(line(pi, q, c)) = inv(line(pi^-1,
q - c, -c)).  Checked exhaustively over all pi, q, c for 4 <= n <= NX_EXHAUST,
and by random sampling for NX_EXHAUST < n <= NX_SAMPLE.

Part B (C38): the exact double-counting identity for the sum, over all n^2
reference points (a, b) in Z_n x Z_n, of the number of "matching" (mutually
non-inverted) pairs of points (i, pi(i)) -- Lemma 4.  Checked the same way
(brute force is O(n^4) per permutation, so the exhaustive range is smaller).

Part C: quantitative record of why the two rejected approaches fail (kept for
the record, not new claims): (i) plain averaging of inv(w) over all n^2 cuts
exceeds floor((n-1)^2/4) already on rev_n; (ii) the single-parameter family
c = 0 (only n candidate cuts) is insufficient for H13-I.

Usage: python3 checks/check_H13I_identities.py [--nx 7] [--ns 40] [--samples 200]
Output: data/runs/check_H13I_identities/report.json, report.md.
"""

import argparse
import itertools
import json
import os
import random
import sys
import time
from fractions import Fraction

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
from moves import CORE_VERSION, sigma, rev  # noqa: E402

VERSION = "check_H13I_identities-1.0"
OUT = os.path.join(ROOT, "data", "runs", "check_H13I_identities")


def line(pi, q, c, n):
    return tuple((pi[(q + 1 + j) % n] - (q + 1 - c)) % n for j in range(n))


def inv(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def inverse_perm(pi):
    n = len(pi)
    out = [0] * n
    for i, v in enumerate(pi):
        out[v] = i
    return tuple(out)


def random_perm(n, rng):
    p = list(range(n))
    rng.shuffle(p)
    return tuple(p)


# ---------------- Part A: shift identities + duality ----------------

def check_c_shift(pi, n):
    for q in range(n):
        for c in range(n):
            w0 = line(pi, q, c, n)
            w1 = line(pi, q, (c + 1) % n, n)
            delta = inv(w1) - inv(w0)
            posmax = w0.index(n - 1)
            pred = 2 * posmax - (n - 1)
            if delta != pred:
                return False, (q, c, delta, pred)
    return True, None


def check_q_shift(pi, n):
    for q in range(n):
        for c in range(n):
            w0 = line(pi, q, c, n)
            w1 = line(pi, (q + 1) % n, c, n)
            delta = inv(w1) - inv(w0)
            v = w0[1:] + w0[:1]
            j0 = v.index(0)
            pred = (n - 1 - 2 * w0[0]) + (n - 1 - 2 * j0)
            if delta != pred:
                return False, (q, c, delta, pred)
    return True, None


def check_duality(pi, n):
    pinv = inverse_perm(pi)
    for q in range(n):
        for c in range(n):
            a = inv(line(pi, q, c, n))
            b = inv(line(pinv, (q - c) % n, (-c) % n, n))
            if a != b:
                return False, (q, c, a, b)
    return True, None


# ---------------- Part B: double counting ----------------

def matches_count(pi, n):
    total = 0
    for a in range(n):
        for b in range(n):
            x = [(i - a) % n for i in range(n)]
            y = [(pi[i] - b) % n for i in range(n)]
            order = sorted(range(n), key=lambda i: x[i])
            m = 0
            for ii in range(n):
                for jj in range(ii + 1, n):
                    if y[order[ii]] < y[order[jj]]:
                        m += 1
            total += m
    return total


def matches_formula(pi, n):
    S = 0
    for i in range(n):
        for ip in range(i + 1, n):
            D = ip - i
            Dp = (pi[ip] - pi[i]) % n
            S += n * n - n * (D + Dp) + 2 * D * Dp
    return S


def check_double_counting(pi, n):
    a = matches_count(pi, n)
    b = matches_formula(pi, n)
    return a == b, (a, b)


# ---------------- Part C: record of rejected approaches ----------------

def avg_inv_all_cuts(pi, n):
    total = 0
    for q in range(n):
        for c in range(n):
            total += inv(line(pi, q, c, n))
    return Fraction(total, n * n)


def part_c(log, c0_exhaustive_max=8):
    """c0_exhaustive_max caps the exhaustive n!-scan for the c=0 family (cheap
    up to n=8, ~4e4 perms; n=9,10 would be 3.6e5/3.6e6 perms x n x O(n^2) each
    -- too slow to redo here and unnecessary, since the family is already
    known insufficient from n <= 8 and from the distinct rejected family in
    docs/notes/h13_line_model.md section 1.7 (n = 5, 7 counterexamples)."""
    rows = []
    for n in range(4, 11):
        bound = (n - 1) ** 2 // 4
        avg_sigma = avg_inv_all_cuts(sigma(n), n)
        avg_rev = avg_inv_all_cuts(rev(n), n)
        c0_worst = 0
        c0_arg = None
        if n <= c0_exhaustive_max:
            for pi in itertools.permutations(range(n)):
                m = min(inv(line(pi, q, 0, n)) for q in range(n))
                if m > c0_worst:
                    c0_worst, c0_arg = m, pi
        else:
            c0_worst, c0_arg = None, "not recomputed (n > %d, see docstring)" % c0_exhaustive_max
        row = {"n": n, "bound_floor((n-1)^2/4)": bound,
               "avg_inv_sigma_n": str(avg_sigma), "avg_inv_rev_n": str(avg_rev),
               "avg_exceeds_bound": bool(avg_rev > bound),
               "max_pi_min_q_inv_c0": c0_worst,
               "c0_family_sufficient": (bool(c0_worst <= bound) if c0_worst is not None else None),
               "argmax_c0_family": (list(c0_arg) if isinstance(c0_arg, tuple) else c0_arg)}
        rows.append(row)
        log(f"n={n}: avg(sigma_n)={float(avg_sigma):.2f} avg(rev_n)={float(avg_rev):.2f} "
            f"bound={bound} -> averaging {'FAILS' if avg_rev > bound else 'holds'}; "
            f"c=0 family worst-case min = {c0_worst}")
    return rows


def run(nx, ns, samples, log):
    t0 = time.time()
    rng = random.Random(12345)
    rows = []
    for n in range(4, nx + 1):
        t1 = time.time()
        fails = {"c_shift": 0, "q_shift": 0, "duality": 0, "double_counting": 0}
        checked = 0
        for pi in itertools.permutations(range(n)):
            checked += 1
            ok, info = check_c_shift(pi, n)
            if not ok:
                fails["c_shift"] += 1
                log(f"FAIL c-shift n={n} pi={pi} info={info}")
            ok, info = check_q_shift(pi, n)
            if not ok:
                fails["q_shift"] += 1
                log(f"FAIL q-shift n={n} pi={pi} info={info}")
            ok, info = check_duality(pi, n)
            if not ok:
                fails["duality"] += 1
                log(f"FAIL duality n={n} pi={pi} info={info}")
            if n <= 6:
                ok, info = check_double_counting(pi, n)
                if not ok:
                    fails["double_counting"] += 1
                    log(f"FAIL double-counting n={n} pi={pi} info={info}")
        row = {"n": n, "mode": "exhaustive", "count": checked, "fails": fails,
               "seconds": round(time.time() - t1, 1)}
        rows.append(row)
        log(f"n={n}: exhaustive, {checked} perms, fails={fails}, {row['seconds']}s")
    for n in range(nx + 1, ns + 1):
        fails = {"c_shift": 0, "q_shift": 0, "duality": 0, "double_counting": 0}
        t1 = time.time()
        for _ in range(samples):
            pi = random_perm(n, rng)
            ok, info = check_c_shift(pi, n)
            if not ok:
                fails["c_shift"] += 1
                log(f"FAIL c-shift n={n} pi={pi} info={info}")
            ok, info = check_q_shift(pi, n)
            if not ok:
                fails["q_shift"] += 1
                log(f"FAIL q-shift n={n} pi={pi} info={info}")
            ok, info = check_duality(pi, n)
            if not ok:
                fails["duality"] += 1
                log(f"FAIL duality n={n} pi={pi} info={info}")
            if n <= 14:
                ok, info = check_double_counting(pi, n)
                if not ok:
                    fails["double_counting"] += 1
                    log(f"FAIL double-counting n={n} pi={pi} info={info}")
        # also check sigma_n and rev_n specifically (structured, not random)
        for name, pi in (("sigma_n", sigma(n)), ("rev_n", rev(n))):
            ok, info = check_c_shift(pi, n)
            if not ok:
                fails["c_shift"] += 1
                log(f"FAIL c-shift n={n} {name} info={info}")
            ok, info = check_q_shift(pi, n)
            if not ok:
                fails["q_shift"] += 1
                log(f"FAIL q-shift n={n} {name} info={info}")
            ok, info = check_duality(pi, n)
            if not ok:
                fails["duality"] += 1
                log(f"FAIL duality n={n} {name} info={info}")
        row = {"n": n, "mode": f"sampled ({samples} random + sigma_n + rev_n)",
               "count": samples + 2, "fails": fails, "seconds": round(time.time() - t1, 1)}
        rows.append(row)
        log(f"n={n}: sampled ({samples} random + sigma_n, rev_n), fails={fails}, {row['seconds']}s")
    total_fails = sum(sum(r["fails"].values()) for r in rows)
    part_c_rows = part_c(log)
    log(f"TOTAL: {'PASS' if total_fails == 0 else 'FAIL'} ({time.time() - t0:.0f} s)")
    return rows, part_c_rows, total_fails


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nx", type=int, default=7, help="exhaustive range 4..nx")
    ap.add_argument("--ns", type=int, default=40, help="sampled range nx+1..ns")
    ap.add_argument("--samples", type=int, default=200)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    logf = open(os.path.join(OUT, "report.log"), "a")

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + "\n")
        logf.flush()

    log(f"== {VERSION} {CORE_VERSION} args={vars(args)}")
    rows, part_c_rows, total_fails = run(args.nx, args.ns, args.samples, log)
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump({"version": VERSION, "core": CORE_VERSION, "args": vars(args),
                    "identity_checks": rows, "rejected_approaches": part_c_rows,
                    "total_fails": total_fails}, f, indent=1)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write(f"# check_H13I_identities report ({VERSION}, {CORE_VERSION})\n\n")
        f.write(f"Overall: {'PASS' if total_fails == 0 else 'FAIL'} (total_fails={total_fails})\n\n")
        f.write("## Part A/B: shift identities (C37), duality (C37), double counting (C38)\n\n")
        f.write("| n | mode | count | c_shift fails | q_shift fails | duality fails | "
                "double_counting fails | seconds |\n|---|---|---|---|---|---|---|---|\n")
        for r in rows:
            f.write(f"| {r['n']} | {r['mode']} | {r['count']} | {r['fails']['c_shift']} | "
                    f"{r['fails']['q_shift']} | {r['fails']['duality']} | "
                    f"{r['fails']['double_counting']} | {r['seconds']} |\n")
        f.write("\n## Part C: why averaging and the c=0 family fail (record, no new claim)\n\n")
        f.write("| n | bound | avg inv sigma_n | avg inv rev_n | averaging fails? | "
                "max_pi min_q inv(.,q,0) | c=0 family sufficient? |\n|---|---|---|---|---|---|---|\n")
        for r in part_c_rows:
            f.write(f"| {r['n']} | {r['bound_floor((n-1)^2/4)']} | {r['avg_inv_sigma_n']} | "
                    f"{r['avg_inv_rev_n']} | {r['avg_exceeds_bound']} | "
                    f"{r['max_pi_min_q_inv_c0']} | {r['c0_family_sufficient']} |\n")
    print(f"PASS" if total_fails == 0 else "FAIL")


if __name__ == "__main__":
    main()
