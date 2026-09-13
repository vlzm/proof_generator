"""Independent numerical checks of H10 / C27 (CLAIMS.md): for the carrier-route
variant of N1 (constructions/strict_upper.py, route="carrier"),

    min_c [2F_c - S_c + R_c] <= B_n + 1   for all n >= 4, all pi in S_n.

This module does not reimplement cycles_of/cycle_data/choose_carrier/
carrier_route from scratch (AGENTS.md rule 3 exemption is for the oracle
move semantics, not this variant, which is itself only measured, not proved
-- rule 16). Instead it uses a *lean* path that skips building and verifying
the actual word (only F_c, S_c, R_c are needed for this inequality), and
cross-validates that lean path against the fully-asserted `shift_word`
(which does build and verify the word letter by letter against oracle-1.0)
on an exhaustive sweep of a smaller n before trusting it on larger n --
see `_selfcheck_against_shift_word`.

Two independent pieces of evidence:

1. `reflections(max_n)`: exhaustive over the reflection family pi(i) = h-i
   mod n (all h, all c), which C24/C25/C27v studied at 4<=n<=12. Finds, for
   each n, the worst h (max over h of min over c of 2F_c-S_c+R_c) and its
   gap over B_n.
2. `perms10()`: exhaustive over all 10! = 3628800 permutations of S_10 (all
   n=10 shifts each), extending the exhaustive range of C27v from n<=9 to
   n=10. Checkpointed (AGENTS.md rule 10): progress is written periodically
   to OUT_DIR/checkpoint.json and the run resumes from there.

Usage:
    python3 checks/check_C27.py --selfcheck            # cross-validate lean path, n<=7
    python3 checks/check_C27.py --reflections 150       # reflection family, 4<=n<=150
    python3 checks/check_C27.py --perms10               # exhaustive n=10 (checkpointed, ~20 min)
Output: data/runs/check_C27/*.json, report.md
"""

import argparse
import itertools
import json
import os
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, ROOT)

from constructions.strict_upper import (  # noqa: E402
    cycles_of, cycle_data, choose_carrier, carrier_route, shift_word, CONSTRUCTION_VERSION,
)
from oracle.moves import CORE_VERSION  # noqa: E402
from bounds.known import target_diameter as B_n  # noqa: E402

OUT_DIR = os.path.join(ROOT, "data", "runs", "check_C27")


def min_expr_lean(pi):
    """min_c [2F_c - S_c + R_c], computed without building/verifying the word."""
    n = len(pi)
    best = None
    for c in range(n):
        f = [(pi[i] + c) % n for i in range(n)]
        F = S = 0
        required = []
        for cyc in cycles_of(f):
            cd = cycle_data(cyc, n)
            F += cd["F"]
            S += cd["S"]
            required.append(cyc[choose_carrier(cyc, cd, n)])
        _, letters = carrier_route(required, c, n)
        expr = 2 * F - S + len(letters)
        if best is None or expr < best:
            best = expr
    return best


def min_expr_asserted(pi):
    """Same quantity via shift_word (builds and asserts the full word). Reference
    for cross-validation, too slow for exhaustive n=10."""
    n = len(pi)
    best = None
    for c in range(n):
        sw = shift_word(list(pi), c, route="carrier")
        st = sw["stats"]
        expr = 2 * st["F_c"] - st["S_c"] + st["route_len"]
        if best is None or expr < best:
            best = expr
    return best


def selfcheck(max_n=7):
    """Exhaustive cross-validation of min_expr_lean against min_expr_asserted."""
    results = {}
    for n in range(4, max_n + 1):
        t0 = time.time()
        checked = mismatches = 0
        for pi in itertools.permutations(range(n)):
            a = min_expr_lean(list(pi))
            b = min_expr_asserted(pi)
            checked += 1
            if a != b:
                mismatches += 1
                if mismatches <= 5:
                    print("MISMATCH", pi, "lean=", a, "asserted=", b)
        results[n] = {"checked": checked, "mismatches": mismatches, "seconds": round(time.time() - t0, 3)}
        print(f"n={n}: {checked} perms, {mismatches} mismatches, {results[n]['seconds']} s")
    return results


def reflections(max_n):
    """Exhaustive over pi(i) = h-i mod n, all h, all c, 4<=n<=max_n."""
    rows = []
    for n in range(4, max_n + 1):
        worst = None
        for h in range(n):
            pi = [(h - i) % n for i in range(n)]
            v = min_expr_lean(pi)
            if worst is None or v > worst[0]:
                worst = (v, h)
        bn = B_n(n)
        rows.append({"n": n, "B_n": bn, "worst_min_expr": worst[0], "diff": worst[0] - bn, "worst_h": worst[1]})
    return rows


def perms10(checkpoint_every=200000):
    """Exhaustive n=10, checkpointed. Returns final summary dict."""
    n = 10
    os.makedirs(OUT_DIR, exist_ok=True)
    ckpt_path = os.path.join(OUT_DIR, "checkpoint_n10.json")
    start_index = 0
    max_diff = None
    worst = []  # list of (diff, pi), largest first, capped
    t_prev_elapsed = 0.0
    if os.path.exists(ckpt_path):
        with open(ckpt_path) as f:
            ck = json.load(f)
        start_index = ck["completed"]
        max_diff = ck["max_diff"]
        worst = ck["worst"]
        t_prev_elapsed = ck.get("elapsed_seconds", 0.0)
        print(f"resuming perms10 from index {start_index} (max_diff so far {max_diff})")
    bn = B_n(n)
    t0 = time.time()
    gen = itertools.permutations(range(n))
    for _ in range(start_index):
        next(gen)
    idx = start_index
    for pi in gen:
        v = min_expr_lean(list(pi))
        diff = v - bn
        if max_diff is None or diff > max_diff:
            max_diff = diff
        if diff >= 0:
            worst.append((diff, list(pi)))
            worst.sort(key=lambda t: -t[0])
            del worst[20:]
        idx += 1
        if idx % checkpoint_every == 0:
            elapsed = t_prev_elapsed + (time.time() - t0)
            with open(ckpt_path, "w") as f:
                json.dump({"completed": idx, "total": 3628800, "max_diff": max_diff,
                           "worst": worst, "elapsed_seconds": round(elapsed, 1)}, f)
            print(f"  checkpoint: {idx}/3628800, max_diff={max_diff}, elapsed={elapsed:.0f}s")
    elapsed = t_prev_elapsed + (time.time() - t0)
    result = {"n": n, "B_n": bn, "count": idx, "max_diff": max_diff, "worst": worst[:20],
              "elapsed_seconds": round(elapsed, 1), "coverage": "all 10! permutations, all 10 shifts, carrier route",
              "construction_version": CONSTRUCTION_VERSION, "core_version": CORE_VERSION}
    with open(os.path.join(OUT_DIR, "result_n10.json"), "w") as f:
        json.dump(result, f, indent=2)
    if os.path.exists(ckpt_path):
        os.remove(ckpt_path)
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selfcheck", type=int, nargs="?", const=7, default=None,
                     help="cross-validate lean path vs shift_word, exhaustive 4<=n<=N (default 7)")
    ap.add_argument("--reflections", type=int, default=0, help="reflection family, 4<=n<=N")
    ap.add_argument("--perms10", action="store_true", help="exhaustive n=10 (checkpointed)")
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)

    if args.selfcheck is not None:
        res = selfcheck(args.selfcheck)
        with open(os.path.join(OUT_DIR, "selfcheck.json"), "w") as f:
            json.dump(res, f, indent=2)

    if args.reflections:
        rows = reflections(args.reflections)
        with open(os.path.join(OUT_DIR, "reflections.json"), "w") as f:
            json.dump(rows, f, indent=2)
        bad = [r for r in rows if r["diff"] not in (0, 1)]
        maxdiff = max(r["diff"] for r in rows)
        print(f"reflections 4<=n<={args.reflections}: max diff over all n = {maxdiff}; "
              f"n with diff==1: {[r['n'] for r in rows if r['diff'] == 1][:10]}...")
        if bad:
            print("UNEXPECTED diff values:", bad)

    if args.perms10:
        res = perms10()
        print(f"perms10: count={res['count']}, max_diff={res['max_diff']}, elapsed={res['elapsed_seconds']}s")


if __name__ == "__main__":
    main()
