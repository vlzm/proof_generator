"""Checker for C37 (session 9): two negative results about H13-I (toric
minimum-inversions bound `I(pi) <= floor((n-1)^2/4)`, docs/notes/h13_line_model.md).

  Part A (exact cut-sum formula). For pi in S_n, define for the double cut
    (q, c) the line w_j = pi(q+1+j) - (q+1-c) mod n, j = 0..n-1 (PROBLEM
    convention, docs/notes/h13_line_model.md §0). Then
      sum_{q,c} inv(w_{q,c}) = sum_{i<i'} [ n(delta+delta') - 2*delta*delta' ]
    where delta = i'-i, delta' = (pi(i')-pi(i)) mod n. Checked against brute
    force over all n^2 cuts, all pi, 4 <= n <= AMAX.

    For reflections pi_h(i) = h-i mod n, delta' = n-delta identically, so the
    formula reduces to a closed form independent of h:
      avg_inv(pi_h) = (1/n^2) sum_{q,c} inv = n(n-1)/2 - (n^2-1)/6
    which is ~n^2/3, a factor 4/3 above the target floor((n-1)^2/4) ~ n^2/4
    for large n (I(pi_h) = floor((n-1)^2/4) exactly, C33). Since min <= avg
    is the only inequality plain averaging over cuts gives, and here
    avg > target by a growing additive amount, uniform averaging over the
    n^2 cuts cannot prove H13-I in general (sharpens the qualitative note
    in h13_line_model.md, session 8, "usredneniye ne rabotayet").

  Part B (single-axis rotation is insufficient). Restricting to only q
  (c = 0 fixed) or only c (q = -1 fixed, i.e. natural start) does not reach
  floor((n-1)^2/4) for some pi starting at n = 7: exhaustive over all pi,
  4 <= n <= AMAX1D, reports max_pi min_q_only(inv) and max_pi min_c_only(inv)
  against the target. Both exceed the target by 1 at n = 7, 8 (example
  (0,5,4,3,2,1,6) at n=7, already the Lemma C counterexample). Confirms the
  double cut genuinely needs both parameters, not a reduction to a known
  1-D cyclic-rotation-inversions lemma.

Neither part proves H13-I; both are negative/structural results ruling out
two natural proof strategies (see docs/notes/h13_line_model.md session 9).

Usage: python3 checks/check_C37.py [--amax 9] [--amax1d 8]
Output: data/runs/check_C37/report.json, report.md.
"""

import argparse
import itertools
import json
import os
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
VERSION = "check_C37-1.0"
OUT = os.path.join(ROOT, "data", "runs", "check_C37")


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def brute_cut_total(pi):
    n = len(pi)
    total = 0
    for q in range(n):
        for c in range(n):
            w = [(pi[(q + 1 + j) % n] - (q + 1 - c)) % n for j in range(n)]
            total += inversions(w)
    return total


def formula_cut_total(pi):
    n = len(pi)
    total = 0
    for i in range(n):
        for ip in range(i + 1, n):
            delta = ip - i
            deltap = (pi[ip] - pi[i]) % n
            total += n * (delta + deltap) - 2 * delta * deltap
    return total


def reflection(n, h):
    return tuple((h - i) % n for i in range(n))


def part_a(amax, amax_reflection, log):
    ok = True
    rows = []
    for n in range(4, amax + 1):
        t0 = time.time()
        checked = 0
        for pi in itertools.permutations(range(n)):
            bt = brute_cut_total(pi)
            ft = formula_cut_total(pi)
            if bt != ft:
                ok = False
                log(f"FAIL part A: n={n} pi={pi} brute={bt} formula={ft}")
            checked += 1
        log(f"part A (all pi): n={n} checked={checked} perms ({time.time()-t0:.1f}s)")

    for n in range(4, amax_reflection + 1):
        t0 = time.time()
        h = n - 1
        pi = reflection(n, h)
        bt = brute_cut_total(pi)
        avg_brute = bt / (n * n)
        avg_formula = n * (n - 1) / 2 - (n * n - 1) / 6
        target = (n - 1) ** 2 // 4
        if abs(avg_brute - avg_formula) > 1e-6:
            ok = False
            log(f"FAIL part A reflection avg: n={n} brute={avg_brute} formula={avg_formula}")
        rows.append(
            {
                "n": n,
                "reflection_avg_brute": avg_brute,
                "reflection_avg_formula": avg_formula,
                "target_floor((n-1)^2/4)": target,
                "avg_over_target": avg_formula / target if target else None,
                "seconds": time.time() - t0,
            }
        )
        log(
            f"part A (reflection avg): n={n} brute={avg_brute:.4f} formula={avg_formula:.4f} "
            f"target={target} ratio={avg_formula/target if target else float('nan'):.4f} ({time.time()-t0:.1f}s)"
        )
    return ok, rows


def min_q_only(pi):
    n = len(pi)
    best = None
    for q in range(n):
        w = [pi[(q + 1 + j) % n] for j in range(n)]
        iv = inversions(w)
        if best is None or iv < best:
            best = iv
    return best


def min_c_only(pi):
    n = len(pi)
    best = None
    for c in range(n):
        w = [(pi[j] + c) % n for j in range(n)]
        iv = inversions(w)
        if best is None or iv < best:
            best = iv
    return best


def part_b(amax1d, log):
    rows = []
    example = None
    for n in range(4, amax1d + 1):
        t0 = time.time()
        target = (n - 1) ** 2 // 4
        worst_q = 0
        worst_c = 0
        for pi in itertools.permutations(range(n)):
            mq = min_q_only(pi)
            mc = min_c_only(pi)
            if mq > worst_q:
                worst_q = mq
            if mc > worst_c:
                worst_c = mc
            if n == 7 and mq > target and example is None:
                example = list(pi)
        rows.append(
            {
                "n": n,
                "target": target,
                "max_min_q_only": worst_q,
                "max_min_c_only": worst_c,
                "seconds": time.time() - t0,
            }
        )
        log(
            f"part B: n={n} target={target} max_min_q_only={worst_q} "
            f"max_min_c_only={worst_c} ({time.time()-t0:.1f}s)"
        )
    return rows, example


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amax", type=int, default=7, help="exhaustive-over-all-pi n for part A (cut-sum formula)")
    ap.add_argument("--amax_reflection", type=int, default=40, help="max n for the reflection closed-form check (single pi per n)")
    ap.add_argument("--amax1d", type=int, default=8, help="exhaustive n for part B (1-D rotation)")
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    lines = []

    def log(msg):
        print(msg)
        lines.append(msg)

    log(f"{VERSION} amax={args.amax} amax_reflection={args.amax_reflection} amax1d={args.amax1d}")
    ok_a, rows_a = part_a(args.amax, args.amax_reflection, log)
    rows_b, example = part_b(args.amax1d, log)

    result = {
        "version": VERSION,
        "part_a_ok": ok_a,
        "part_a_rows": rows_a,
        "part_b_rows": rows_b,
        "part_b_example_n7": example,
    }
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(result, f, indent=2)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# check_C37 report\n\n" + "\n".join(lines) + "\n")
        f.write(f"\nPart A overall: {'PASS' if ok_a else 'FAIL'}\n")
        f.write(f"Part B example (n=7, min_q_only > target): {example}\n")

    print("PASS" if ok_a else "FAIL")


if __name__ == "__main__":
    main()
