"""Checker for the H13-I session (session 9): the proved lemmas of
docs/notes/h13i_toric.md and the exhaustive data of data/runs/h13i_toric/.

I(pi) = min over all n^2 double cuts (a, b) of inv(w), w_j = pi(a + j) - b (mod n)
      (a = q + 1, b = q + 1 - c in the (q, c) notation of h13_line_model.md).

Parts (each independent of experiments/toric_inv.c: every quantity here is
recomputed from the definition by brute force):

  A. Lemma 1 (shift recurrences).  For all pi and all cuts, 4 <= n <= AMAX:
       inv(a+1, b) = inv(a, b) + n - 1 - 2*((pi(a) - b) mod n)
       inv(a, b+1) = inv(a, b) + n - 1 - 2*((pi^{-1}(b) - a) mod n)
  B. Lemma 3 (reversal duality).  I(pi) + M(pi^r) = C(n, 2) for all pi, where
     pi^r(i) = pi(-i) and M = max over the n^2 cuts; plus the arithmetic
     identity floor(n^2/4) + floor((n-1)^2/4) = C(n, 2) for 4 <= n <= 200.
     Hence H13-I is equivalent to "every pi has a cut with >= floor(n^2/4)
     inversions", with the equality cases swapped (reflections <-> rotations).
  C. Lemma 4 (winding number).  omega(pi) = (1/n) sum_i ((pi(i+1)-pi(i)) mod n)
     is an integer in [1, n-1], invariant under both rotations; for EVERY value
     cut b the cyclic line has exactly omega descents; omega = 1 iff pi is a
     rotation; omega = n-1 iff pi is a reflection.
  D. Lemma 5 (reflections).  I(pi_h) = floor((n-1)^2/4) for every reflection
     pi_h(i) = (h - i) mod n: brute force for 4 <= n <= 9, and the closed form
     min_k [C(k,2) + C(n-k,2)] = floor((n-1)^2/4) for 4 <= n <= 200.
  E. Lemma 6 (increasing runs).  Cutting at a cyclic descent gives a line that
     is a concatenation of exactly omega increasing runs, whence
     inv <= (n^2 - sum n_L^2)/2; in particular I = 0 when omega = 1 and
     I <= floor(n^2/4) when omega = 2.
  G. Lemma 2 (pair areas).  A pair of points at cyclic distances (dx, dy) is
     concordant at exactly dx*dy + (n-dx)*(n-dy) of the n^2 cuts, so the total
     concordance over all cuts is n * sum_delta delta * W_delta, where
     W_delta = #{i : pi(i+delta) < pi(i)} and sum_delta W_delta = C(n, 2).
  F. H13-I itself.  Exhaustive brute force I(pi) <= floor((n-1)^2/4) for
     4 <= n <= AMAX with equality exactly on the n reflections; cross-check of
     data/runs/h13i_toric/toric_n*.json (toric_inv-1.0, 4 <= n <= 13) and of the
     independent session-8 histograms data/runs/line_profile/profile_n*.json
     (line_profile-1.0, 4 <= n <= 10): full I-histograms must agree.

Usage: python3 checks/check_H13I.py [--amax 7]
Output: data/runs/check_H13I/report.json, report.md.  Exit code 0 on PASS.
"""

import argparse
import itertools
import json
import os
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

VERSION = "check_H13I-1.0"
OUT = os.path.join(ROOT, "data", "runs", "check_H13I")


def target(n):
    return ((n - 1) ** 2) // 4


def line(pi, n, a, b):
    return [(pi[(a + j) % n] - b) % n for j in range(n)]


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def cut_table(pi, n):
    """All n^2 inversion counts, straight from the definition."""
    return [[inversions(line(pi, n, a, b)) for b in range(n)] for a in range(n)]


def winding(pi, n):
    return sum((pi[(i + 1) % n] - pi[i]) % n for i in range(n)) // n


def is_rotation(pi, n):
    return all((pi[i] - i) % n == (pi[0] - 0) % n for i in range(n))


def is_reflection(pi, n):
    return all((pi[i] + i) % n == pi[0] % n for i in range(n))


def part_a(amax, log):
    ok = True
    checked = 0
    for n in range(4, amax + 1):
        for pi in itertools.permutations(range(n)):
            t = cut_table(pi, n)
            ipi = [0] * n
            for i, v in enumerate(pi):
                ipi[v] = i
            for a in range(n):
                for b in range(n):
                    lhs = t[(a + 1) % n][b]
                    rhs = t[a][b] + n - 1 - 2 * ((pi[a] - b) % n)
                    if lhs != rhs:
                        ok = False
                        log.append(f"A FAIL n={n} pi={pi} a={a} b={b} (position shift)")
                    lhs = t[a][(b + 1) % n]
                    rhs = t[a][b] + n - 1 - 2 * ((ipi[b] - a) % n)
                    if lhs != rhs:
                        ok = False
                        log.append(f"A FAIL n={n} pi={pi} a={a} b={b} (value shift)")
                    checked += 2
    log.append(f"A: shift recurrences verified on {checked} cut pairs, 4 <= n <= {amax}: "
               + ("PASS" if ok else "FAIL"))
    return ok, {"checked_transitions": checked, "amax": amax}


def part_b(amax, log):
    ok = True
    cnt = 0
    for n in range(4, amax + 1):
        for pi in itertools.permutations(range(n)):
            rev = [pi[(-i) % n] for i in range(n)]
            t = cut_table(pi, n)
            tr = cut_table(rev, n)
            I = min(min(r) for r in t)
            M = max(max(r) for r in tr)
            if I + M != n * (n - 1) // 2:
                ok = False
                log.append(f"B FAIL n={n} pi={pi}: I={I} M(pi^r)={M}")
            cnt += 1
    ident = all(((n * n) // 4) + target(n) == n * (n - 1) // 2 for n in range(4, 201))
    if not ident:
        ok = False
        log.append("B FAIL: floor(n^2/4) + floor((n-1)^2/4) != C(n,2) somewhere in 4..200")
    log.append(f"B: duality I(pi) + M(pi^r) = C(n,2) on {cnt} permutations (4 <= n <= {amax}); "
               f"arithmetic identity 4 <= n <= 200: " + ("PASS" if ok else "FAIL"))
    return ok, {"permutations": cnt, "amax": amax, "arith_identity_upto": 200}


def part_c(amax, log):
    ok = True
    cnt = 0
    for n in range(4, amax + 1):
        for pi in itertools.permutations(range(n)):
            om = winding(pi, n)
            if not (1 <= om <= n - 1):
                ok = False
                log.append(f"C FAIL n={n} pi={pi}: omega={om} out of range")
            # invariance under both rotations
            for s in range(n):
                rot_pos = [pi[(i + s) % n] for i in range(n)]
                rot_val = [(pi[i] + s) % n for i in range(n)]
                if winding(rot_pos, n) != om or winding(rot_val, n) != om:
                    ok = False
                    log.append(f"C FAIL n={n} pi={pi}: omega not rotation-invariant")
            # exactly omega cyclic descents for every value cut b
            for b in range(n):
                d = sum(1 for i in range(n)
                        if (pi[(i + 1) % n] - b) % n < (pi[i] - b) % n)
                if d != om:
                    ok = False
                    log.append(f"C FAIL n={n} pi={pi} b={b}: {d} cyclic descents, omega={om}")
            if (om == 1) != is_rotation(pi, n):
                ok = False
                log.append(f"C FAIL n={n} pi={pi}: omega=1 vs rotation mismatch")
            if (om == n - 1) != is_reflection(pi, n):
                ok = False
                log.append(f"C FAIL n={n} pi={pi}: omega=n-1 vs reflection mismatch")
            cnt += 1
    log.append(f"C: winding lemma on {cnt} permutations, 4 <= n <= {amax}: "
               + ("PASS" if ok else "FAIL"))
    return ok, {"permutations": cnt, "amax": amax}


def part_d(log):
    ok = True
    for n in range(4, 10):
        for h in range(n):
            pi = [(h - i) % n for i in range(n)]
            t = cut_table(pi, n)
            I = min(min(r) for r in t)
            if I != target(n):
                ok = False
                log.append(f"D FAIL n={n} h={h}: I={I} != {target(n)}")
    for n in range(4, 201):
        m = min((k * (k - 1)) // 2 + ((n - k) * (n - k - 1)) // 2 for k in range(n + 1))
        if m != target(n):
            ok = False
            log.append(f"D FAIL closed form n={n}: min_k = {m} != {target(n)}")
    log.append("D: reflections I = floor((n-1)^2/4) (brute force 4 <= n <= 9, closed form "
               "4 <= n <= 200): " + ("PASS" if ok else "FAIL"))
    return ok, {"bruteforce_upto": 9, "closed_form_upto": 200}


def part_e(amax, log):
    ok = True
    cnt = 0
    for n in range(4, amax + 1):
        for pi in itertools.permutations(range(n)):
            om = winding(pi, n)
            for b in range(n):
                for a in range(n):
                    # cut at a cyclic descent: the edge (a-1, a) must be a descent
                    prev = (pi[(a - 1) % n] - b) % n
                    cur = (pi[a] - b) % n
                    if not (cur < prev):
                        continue
                    w = line(pi, n, a, b)
                    runs = []
                    cur_len = 1
                    for j in range(1, n):
                        if w[j] > w[j - 1]:
                            cur_len += 1
                        else:
                            runs.append(cur_len)
                            cur_len = 1
                    runs.append(cur_len)
                    if len(runs) != om:
                        ok = False
                        log.append(f"E FAIL n={n} pi={pi} a={a} b={b}: {len(runs)} runs, omega={om}")
                    bound = (n * n - sum(r * r for r in runs)) // 2
                    if inversions(w) > bound:
                        ok = False
                        log.append(f"E FAIL n={n} pi={pi} a={a} b={b}: inv > run bound")
                    cnt += 1
    log.append(f"E: run decomposition at a descent on {cnt} (pi, cut) pairs, 4 <= n <= {amax}: "
               + ("PASS" if ok else "FAIL"))
    return ok, {"pairs": cnt, "amax": amax}


def part_f(amax, log):
    ok = True
    hist_ours = {}
    for n in range(4, amax + 1):
        h = {}
        eq = []
        for pi in itertools.permutations(range(n)):
            I = min(min(r) for r in cut_table(pi, n))
            h[I] = h.get(I, 0) + 1
            if I > target(n):
                ok = False
                log.append(f"F FAIL n={n} pi={pi}: I={I} > {target(n)}")
            elif I == target(n):
                eq.append(pi)
        hist_ours[n] = h
        if len(eq) != n or not all(is_reflection(p, n) for p in eq):
            ok = False
            log.append(f"F FAIL n={n}: equality set is not exactly the n reflections ({len(eq)})")
    log.append(f"F1: brute force H13-I with equality exactly on reflections, 4 <= n <= {amax}: "
               + ("PASS" if ok else "FAIL"))

    # cross-check with experiments/toric_inv.c output
    tor = {}
    d = os.path.join(ROOT, "data", "runs", "h13i_toric")
    for n in range(4, 14):
        p = os.path.join(d, f"toric_n{n}.json")
        if not os.path.exists(p):
            continue
        j = json.load(open(p))
        tor[n] = j
        if j["max_I"] != target(n) or not j["holds"] or j["count_I_gt_target"] != 0:
            ok = False
            log.append(f"F FAIL toric_inv n={n}: max_I={j['max_I']} target={target(n)}")
        if j["count_I_eq_max_all"] != n or not j["argmax_all_reflections"]:
            ok = False
            log.append(f"F FAIL toric_inv n={n}: equality set size {j['count_I_eq_max_all']}")
        if n in hist_ours:
            hh = {int(k): v * n for k, v in j["hist_I_normalised"].items()}
            if hh != hist_ours[n]:
                ok = False
                log.append(f"F FAIL toric_inv n={n}: histogram differs from brute force")
    log.append(f"F2: toric_inv-1.0 JSONs for n in {sorted(tor)}: "
               + ("PASS" if ok else "FAIL"))

    # cross-check with the independent session-8 implementation (line_profile-1.0)
    lp = []
    for n in range(4, 11):
        p = os.path.join(ROOT, "data", "runs", "line_profile", f"profile_n{n}.json")
        if not os.path.exists(p):
            continue
        j = json.load(open(p))
        hist8 = {r["I"]: r["count"] for r in j["by_I"]}
        if n in tor:
            hh = {int(k): v * n for k, v in tor[n]["hist_I_normalised"].items()}
            if hh != hist8:
                ok = False
                log.append(f"F FAIL line_profile n={n}: histogram differs from toric_inv")
        if max(hist8) != target(n):
            ok = False
            log.append(f"F FAIL line_profile n={n}: max I = {max(hist8)} != {target(n)}")
        lp.append(n)
    log.append(f"F3: session-8 line_profile-1.0 histograms agree for n in {lp}: "
               + ("PASS" if ok else "FAIL"))
    return ok, {"bruteforce_amax": amax, "toric_inv_n": sorted(tor), "line_profile_n": lp}


def part_g(amax, log):
    """Lemma 2: a pair {p, p'} is concordant at exactly dx*dy + (n-dx)*(n-dy) of the
    n^2 cuts; hence sum over cuts of #concordant = n * sum_delta delta * W_delta,
    and sum_delta W_delta = C(n, 2) for every pi."""
    ok = True
    cnt = 0
    for n in range(4, amax + 1):
        for pi in itertools.permutations(range(n)):
            t = cut_table(pi, n)
            total = sum(n * (n - 1) // 2 - t[a][b] for a in range(n) for b in range(n))
            pair = 0
            for i in range(n):
                for j in range(i + 1, n):
                    dx = (i - j) % n
                    dy = (pi[i] - pi[j]) % n
                    pair += dx * dy + (n - dx) * (n - dy)
            W = [sum(1 for i in range(n) if pi[(i + d) % n] < pi[i]) for d in range(n)]
            wsum = n * sum(d * W[d] for d in range(1, n))
            if not (total == pair == wsum):
                ok = False
                log.append(f"G FAIL n={n} pi={pi}: {total} {pair} {wsum}")
            if sum(W[1:]) != n * (n - 1) // 2:
                ok = False
                log.append(f"G FAIL n={n} pi={pi}: sum W_delta = {sum(W[1:])}")
            cnt += 1
    log.append(f"G: pair-area and average-concordance formulas on {cnt} permutations, "
               f"4 <= n <= {amax}: " + ("PASS" if ok else "FAIL"))
    return ok, {"permutations": cnt, "amax": amax}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amax", type=int, default=7,
                    help="exhaustive brute-force range for parts A, B, C, E, F1 (default 7)")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    log = []
    t0 = time.time()
    res = {}
    allok = True
    for name, fn in (("A", lambda: part_a(args.amax, log)),
                     ("B", lambda: part_b(args.amax, log)),
                     ("C", lambda: part_c(args.amax, log)),
                     ("D", lambda: part_d(log)),
                     ("E", lambda: part_e(min(args.amax, 6), log)),
                     ("F", lambda: part_f(args.amax, log)),
                     ("G", lambda: part_g(min(args.amax, 6), log))):
        ok, info = fn()
        res[name] = {"pass": ok, **info}
        allok = allok and ok
    secs = time.time() - t0
    report = {"version": VERSION, "amax": args.amax, "seconds": round(secs, 1),
              "result": "PASS" if allok else "FAIL", "parts": res, "log": log}
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(report, f, indent=1)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write(f"# check_H13I — {report['result']}\n\n")
        f.write(f"Версия {VERSION}; brute-force диапазон 4 <= n <= {args.amax}; "
                f"время {secs:.1f} с.\n\n")
        for ln in log:
            f.write(f"- {ln}\n")
    for ln in log:
        print(ln)
    print(f"{VERSION}: {report['result']} ({secs:.1f} s)")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
