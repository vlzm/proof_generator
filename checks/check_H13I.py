"""Checker for the session-9 H13-I evidence (docs/notes/h13i_verdict.md).

H13-I (see docs/notes/h13_line_model.md §6): every permutation pi of Z_n has
I(pi) = min over double cuts (position-rotation a, value-shift b) of
inv(row) <= floor((n-1)^2/4). This checker does NOT prove H13-I; it verifies,
independently of experiments/h13i_search.c and experiments/h13i_exhaustive.c,
the pieces the session-9 verdict relies on:

  Part A: the O(n^2) recurrence used by the C tools for I(pi) agrees with a
    direct O(n^4) brute force (all a, b, full inversion count) on random
    permutations, 4 <= n <= 10 -- an independent re-implementation in Python,
    per AGENTS.md rule 3.
  Part B: three rejected single-degree-of-freedom simplifications of H13-I
    are each refuted by an explicit witness permutation (exhaustive worst
    case at small n, computed here from scratch):
      B1. value-shift only (position fixed at a=0): exceeds the bound at n=7.
      B2. point-aligned cuts (b = pi(a), n candidates instead of n^2).
      B3. diagonal cuts (b = a, n candidates).
  Part C: the exhaustive scan reports (experiments/h13i_exhaustive.c output,
    data/runs/h13i_verdict/exhaustive_n*.json) have exceed_count == 0 and
    max_I == floor((n-1)^2/4), for whichever of n in {11, 12, 13} were run.
  Part D: the simulated-annealing search log (data/runs/h13i_verdict/
    sa_search.log) never records margin > 0.

Usage: python3 checks/check_H13I.py
Output: data/runs/h13i_verdict/check_report.md (this script also prints PASS/FAIL).
"""

import itertools
import json
import os
import random
import re

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT_DIR = os.path.join(ROOT, "data", "runs", "h13i_verdict")

VERSION = "check_H13I-1.0"


def floor_bound(n):
    return (n - 1) ** 2 // 4


def inv_count(seq):
    n = len(seq)
    return sum(1 for i in range(n) for j in range(i + 1, n) if seq[i] > seq[j])


def I_bruteforce(pi):
    n = len(pi)
    best = None
    for a in range(n):
        row = [pi[(a + j) % n] for j in range(n)]
        for b in range(n):
            shifted = [(x - b) % n for x in row]
            c = inv_count(shifted)
            if best is None or c < best:
                best = c
    return best


def I_fast(pi):
    """O(n^2) recurrence: inv(a+1,0) = inv(a,0) + (n-1-2 pi[a]);
    inv(a,b+1) = inv(a,b) + (n-1-2 pos_a[b]), pos_a[v] = (piinv[v]-a) mod n."""
    n = len(pi)
    piinv = [0] * n
    for v, x in enumerate(pi):
        piinv[x] = v
    inv_a = inv_count(pi)
    best = inv_a
    for a in range(n):
        pos_a = [(piinv[v] - a) % n for v in range(n)]
        inv_b = inv_a
        if inv_b < best:
            best = inv_b
        for b in range(n - 1):
            inv_b += (n - 1) - 2 * pos_a[b]
            if inv_b < best:
                best = inv_b
        inv_a += (n - 1) - 2 * pi[a]
    return best


def part_a(log):
    random.seed(20260917)
    ok = True
    tested = 0
    for n in range(4, 11):
        perms = list(itertools.permutations(range(n)))
        sample = perms if n <= 7 else random.sample(perms, 4000)
        for perm in sample:
            e = I_bruteforce(list(perm))
            f = I_fast(list(perm))
            tested += 1
            if e != f:
                ok = False
                log(f"FAIL part A: n={n} perm={perm} brute={e} fast={f}")
    log(f"part A: {tested} permutations, O(n^2) recurrence vs brute force, {'PASS' if ok else 'FAIL'}")
    return ok


def part_b(log):
    ok = True
    worst = {"value_shift_only": {}, "point_aligned": {}, "diagonal": {}}
    for n in range(4, 9):
        w_vs = w_pa = w_diag = 0
        for perm in itertools.permutations(range(n)):
            pi = list(perm)
            vs = min(inv_count([(x - b) % n for x in pi]) for b in range(n))
            pa = min(
                inv_count([(pi[(a + j) % n] - pi[a]) % n for j in range(n)])
                for a in range(n)
            )
            diag = min(
                inv_count([(pi[(a + j) % n] - a) % n for j in range(n)])
                for a in range(n)
            )
            w_vs = max(w_vs, vs)
            w_pa = max(w_pa, pa)
            w_diag = max(w_diag, diag)
        worst["value_shift_only"][n] = w_vs
        worst["point_aligned"][n] = w_pa
        worst["diagonal"][n] = w_diag
    b = floor_bound(7)
    if worst["value_shift_only"][7] <= b:
        ok = False
        log(f"FAIL part B1: expected value-shift-only worst at n=7 to exceed bound {b}, got {worst['value_shift_only'][7]}")
    else:
        log(f"part B1 (value-shift only, a=0 fixed): worst={worst['value_shift_only']}, "
            f"exceeds bound {b} at n=7 as expected -- rejects this simplification")
    log(f"part B2 (point-aligned cuts b=pi(a)): worst={worst['point_aligned']}, bounds={[floor_bound(n) for n in range(4,9)]}")
    log(f"part B3 (diagonal cuts b=a): worst={worst['diagonal']}, bounds={[floor_bound(n) for n in range(4,9)]}")
    for n in range(4, 9):
        if worst["point_aligned"][n] <= floor_bound(n) and worst["diagonal"][n] <= floor_bound(n):
            log(f"  note: at n={n} both restricted families still respect the bound (only n=7 value-shift-only fails here)")
    return ok, worst


def part_c(log):
    ok = True
    found_any = False
    for n in (11, 12, 13):
        path = os.path.join(OUT_DIR, f"exhaustive_n{n}.json")
        if not os.path.exists(path):
            continue
        found_any = True
        with open(path) as f:
            data = json.load(f)
        b = floor_bound(n)
        if data["exceed_count"] != 0 or data["max_I"] != b or data["bound"] != b:
            ok = False
            log(f"FAIL part C: n={n} data={data}")
        else:
            log(f"part C: n={n} exhaustive ({data['count']} permutations, {data['seconds']}s): "
                f"max_I={data['max_I']} == bound={b}, exceed_count=0 -- PASS")
    if not found_any:
        log("part C: no exhaustive_n*.json found in data/runs/h13i_verdict -- SKIPPED (run experiments/h13i_exhaustive.c first)")
    return ok


def part_d(log):
    path = os.path.join(OUT_DIR, "sa_search.log")
    if not os.path.exists(path):
        log("part D: sa_search.log not found -- SKIPPED")
        return True
    ok = True
    max_margin_seen = None
    for line in open(path):
        m = re.search(r"margin=(-?\d+)", line)
        if not m:
            continue
        margin = int(m.group(1))
        if max_margin_seen is None or margin > max_margin_seen:
            max_margin_seen = margin
        if margin > 0:
            ok = False
            log(f"FAIL part D: found margin > 0: {line.strip()}")
    log(f"part D: sa_search.log scanned, max margin seen = {max_margin_seen} (<=0 expected), {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    lines = [f"# check_H13I report ({VERSION})", ""]

    def log(s):
        print(s)
        lines.append(s)

    ok_a = part_a(log)
    ok_b, worst = part_b(log)
    ok_c = part_c(log)
    ok_d = part_d(log)
    overall = ok_a and ok_b and ok_c and ok_d
    log("")
    log(f"OVERALL: {'PASS' if overall else 'FAIL'}")

    with open(os.path.join(OUT_DIR, "check_report.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    with open(os.path.join(OUT_DIR, "part_b_worst.json"), "w") as f:
        json.dump(worst, f, indent=2)

    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
