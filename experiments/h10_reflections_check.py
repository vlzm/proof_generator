"""H10 (C27) on reflections, exhaustive over h and c: extends C27v's verified
range for the family of reflections pi_h(i) = (h-i) mod n well beyond the
n<=12 used in experiments/loss_map.py --families.

Version h10_reflections_check-1.0. Core oracle-1.0, construction
strict_upper-1.0 (route="carrier").

Goal: for every n in [n_min, n_max], every h in [0,n), every c in [0,n),
compute 2*F_c - S_c + R_c (route="carrier") and check
  max_h min_c [2F_c - S_c + R_c] - B_n  in {0, 1},
  == 1  iff n % 4 == 2  (docs/proofs/h10_reflections.md section 4).
No permutation/distance table is needed (R_c, F_c, S_c only; d(pi) not used).

Usage:
  python3 experiments/h10_reflections_check.py --max-n 120
Output: data/runs/strict_loss_audit/h10_reflections_check.json
"""

import argparse
import json
import os
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
for sub in ("oracle", "exact", "bounds", "constructions"):
    sys.path.insert(0, os.path.join(ROOT, sub))

from moves import CORE_VERSION  # noqa: E402
import known  # noqa: E402
import strict_upper as su  # noqa: E402

EXPERIMENT_VERSION = "h10_reflections_check-1.0"
OUT = os.path.join(ROOT, "data", "runs", "strict_loss_audit", "h10_reflections_check.json")


def reflection(h, n):
    return tuple((h - i) % n for i in range(n))


def best_over_c(pi, n):
    best = None
    for c in range(n):
        sw = su.shift_word(pi, c, "carrier")
        st = sw["stats"]
        val = 2 * st["F_c"] - st["S_c"] + st["route_len"]
        if best is None or val < best:
            best = val
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=120)
    args = ap.parse_args()
    t0 = time.time()
    rows = {}
    mismatches = []

    def dump():
        out = {
            "meta": {"experiment": EXPERIMENT_VERSION, "construction": su.CONSTRUCTION_VERSION,
                     "core": CORE_VERSION, "command": " ".join(sys.argv),
                     "range": [args.min_n, args.max_n], "seconds_total": round(time.time() - t0, 1)},
            "claim": "max_h min_c [2F_c - S_c + R_c] - B_n == (1 if n%4==2 else 0), route=carrier",
            "mismatches": mismatches,
            "rows": rows,
        }
        os.makedirs(os.path.dirname(OUT), exist_ok=True)
        with open(OUT, "w") as f:
            json.dump(out, f, indent=1)

    for n in range(args.min_n, args.max_n + 1):
        tn0 = time.time()
        Bn = known.target_diameter(n)
        worst = None
        worst_h = None
        for h in range(n):
            v = best_over_c(reflection(h, n), n)
            if worst is None or v > worst:
                worst, worst_h = v, h
        excess = worst - Bn
        expected = 1 if n % 4 == 2 else 0
        ok = excess == expected
        if not ok:
            mismatches.append(n)
        rows[n] = {"B_n": Bn, "max_h_min_c": worst, "excess": excess,
                   "expected_excess": expected, "match": ok,
                   "argmax_h": worst_h, "seconds": round(time.time() - tn0, 2)}
        print(f"n={n}: excess={excess} expected={expected} {'OK' if ok else 'MISMATCH'} "
              f"({rows[n]['seconds']}s)", flush=True)
        dump()
    print(f"done: {len(mismatches)} mismatches over n={args.min_n}..{args.max_n}, "
          f"{round(time.time() - t0, 1)}s -> {OUT}")


if __name__ == "__main__":
    main()
