"""H10 (C27/PLAN §8) on reflections: exhaustive check and closed-form data for
`min_c [2F_c - S_c + R_c]` where R_c is the carrier-route variant of N1 §3
(`constructions/strict_upper.py: carrier_route`).

Version reflection_route-1.0. Core oracle-1.0; construction strict_upper-1.0.

For a reflection pi(i) = h - i mod n, f_c(i) = (pi(i)+c) mod n = (k-i) mod n
with k = (h+c) mod n, i.e. f_c is itself a reflection with parameter k. This
gives closed forms (self-contained derivation, cross-checked below against
the general, audited machinery of strict_upper.py on 4<=n<=15):

- Fixed points of f_c: solutions of 2i = k (mod n). n odd: exactly one,
  i0 = k * inv2 mod n. n even: two if k even, none if k odd.
- Every non-fixed i pairs with j = (k-i) mod n into a transposition. Let
  m = (j-i) mod n (taking i<j as integers in 0..n-1) and L = min(m, n-m)
  (the shorter-arc length, N1 §1 F-contribution 2L). The carrier (N1 §2.1,
  the vertex where the incoming step is negative and the outgoing positive)
  is i if m <= n//2, else j: it is always the start, in the positive
  direction, of the shorter arc between i and j. Every transposition
  (antipodal or not) has S(C) = 3 (PROBLEM §5.3), so S_c = 3 * (number of
  nontrivial pairs).

R_c is then `carrier_route(carriers, c, n)` unchanged (visits the carrier
set plus 0 and c by the shortest walk on the cycle).

This module does NOT reimplement N1's local words or claim a proof; it only
extends the exhaustive family check of C27v (`experiments/loss_map.py`,
4<=n<=12) to much larger n, since R_c and F_c, S_c need no distance table
and no per-cycle word construction here.
"""

import argparse
import json
import os
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
sys.path.insert(0, os.path.join(ROOT, "constructions"))
sys.path.insert(0, os.path.join(ROOT, "bounds"))

from strict_upper import carrier_route, shift_word, CONSTRUCTION_VERSION  # noqa: E402
from moves import CORE_VERSION  # noqa: E402
import known  # noqa: E402

VERSION = "reflection_route-1.0"
OUT_DIR = os.path.join(ROOT, "data", "runs", "h10_reflection_route")


def reflection_FSc(n, k):
    """F_c, S_c, carrier list of f_c for the reflection pi(i) = h-i (any h with
    h+c = k mod n): see module docstring."""
    carriers = []
    F = 0
    for i in range(n):
        j = (k - i) % n
        if j == i or j < i:
            continue
        m = j - i
        L = min(m, n - m)
        F += 2 * L
        carriers.append(i if m <= n // 2 else j)
    S = 3 * len(carriers)
    return F, S, carriers


def Qmin_R_for_h(n, h):
    """min over c of 2F_c - S_c + R_c for the reflection with parameter h."""
    best = None
    for c in range(n):
        k = (h + c) % n
        F, S, carriers = reflection_FSc(n, k)
        _, letters = carrier_route(carriers, c, n)
        val = 2 * F - S + len(letters)
        if best is None or val < best[0]:
            best = (val, c)
    return best


def self_test(nmax=15):
    """Cross-check reflection_FSc + carrier_route against the general, audited
    strict_upper.shift_word(..., route="carrier") on every (n, h, c)."""
    cnt = 0
    for n in range(4, nmax + 1):
        for h in range(n):
            pi = tuple((h - i) % n for i in range(n))
            for c in range(n):
                k = (h + c) % n
                F, S, carriers = reflection_FSc(n, k)
                sw = shift_word(pi, c, "carrier")
                st = sw["stats"]
                assert F == st["F_c"] and S == st["S_c"], (n, h, c, F, S, st)
                assert sorted(carriers) == sorted(st["carriers"]), (n, h, c)
                cnt += 1
    print(f"self_test OK: {cnt} (n,h,c) triples, 4<=n<={nmax}, "
          f"{CORE_VERSION} {CONSTRUCTION_VERSION} {VERSION}")


def scan(nmin, nmax, log=print):
    """Exhaustive over all reflections (all h) and all shifts c, 4<=n<=nmax.
    Returns {n: {"B_n":..., "worst_excess":..., "worst_h":..., "worst_c":...,
    "excess_by_h": [...]}}."""
    out = {}
    for n in range(nmin, nmax + 1):
        t0 = time.time()
        bn = known.target_diameter(n)
        excess_by_h = []
        worst = None
        for h in range(n):
            val, c = Qmin_R_for_h(n, h)
            exc = val - bn
            excess_by_h.append(exc)
            if worst is None or exc > worst[0]:
                worst = (exc, h, c, val)
        out[n] = {"B_n": bn, "worst_excess": worst[0], "worst_h": worst[1],
                   "worst_c": worst[2], "seconds": round(time.time() - t0, 2),
                   "excess_by_h": excess_by_h}
        log(f"n={n}: B_n={bn}, worst min_c[2F_c-S_c+R_c]-B_n = {worst[0]:+d} "
            f"at h={worst[1]} c={worst[2]} ({out[n]['seconds']}s)")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", type=int, default=15)
    ap.add_argument("--nmin", type=int, default=4)
    ap.add_argument("--nmax", type=int, default=140)
    args = ap.parse_args()
    self_test(args.selftest)
    os.makedirs(OUT_DIR, exist_ok=True)
    t0 = time.time()
    result = scan(args.nmin, args.nmax)
    report = {
        "version": VERSION, "core": CORE_VERSION, "construction": CONSTRUCTION_VERSION,
        "goal": "PLAN §8 (H10/C27): extend the exhaustive reflection check of C27v "
                "(4<=n<=12, experiments/loss_map.py) to larger n using the closed "
                "form of this module (cross-checked against strict_upper.py on "
                f"4<=n<={args.selftest})",
        "coverage": f"all n reflections and all n shifts c, exhaustive, {args.nmin}<=n<={args.nmax}",
        "command": f"python3 experiments/reflection_route.py --nmin {args.nmin} --nmax {args.nmax}",
        "seconds_total": round(time.time() - t0, 1),
        "result": result,
    }
    with open(os.path.join(OUT_DIR, "report.json"), "w") as f:
        json.dump(report, f, indent=1)
    maxexc = max(v["worst_excess"] for v in result.values())
    print(f"max excess over {args.nmin}<=n<={args.nmax}: {maxexc:+d}; report written to {OUT_DIR}/report.json")


if __name__ == "__main__":
    main()
