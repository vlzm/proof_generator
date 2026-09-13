"""Affine family against H10 / C27 (carrier-route variant of construction N1).

Version affine_C27-1.0.  Core oracle-1.0, construction strict_upper-1.0.

For pi(i) = a*i + b mod n with gcd(a, n) = 1 (a = n-1 gives the reflections)
this builds, for every shift c, the word of the carrier-route variant and
records

    len_R  = 2F_c - S_c + R_c = |word|      (the quantity of H10/C27)
    red_R  = |freely reduced word|          (the second half of C27)
    len_1  = 2F_c - S_c + H(c)              (N1's own route, C26/H9)
    red_1  = |freely reduced N1 word|

and reports, per n, the maximum over the family of the minimum over c, minus
B_n = n(n-1)/2.  Every word is executed by the reference moves inside
strict_upper.shift_word, so a reported length is the length of a word that
really sorts pi.

Usage:  python3 experiments/affine_C27.py --max-n 30 [--min-n 4]
Output: data/runs/carrier_route_H10/affine.json  (+ progress on stdout)
"""

import argparse
import json
import os
import sys
import time
from math import gcd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
for sub in ("oracle", "bounds", "constructions"):
    sys.path.insert(0, os.path.join(ROOT, sub))

from moves import freely_reduce, CORE_VERSION  # noqa: E402
import known  # noqa: E402
import strict_upper as su  # noqa: E402

EXPERIMENT_VERSION = "affine_C27-1.0"
OUT_DIR = os.path.join(ROOT, "data", "runs", "carrier_route_H10")


def analyze(pi):
    """min over c of the four lengths, and the argmin shifts."""
    n = len(pi)
    best = {}
    for c in range(n):
        wR = su.shift_word(pi, c, "carrier")["word"]
        w1 = su.shift_word(pi, c, "n1")["word"]
        vals = {"len_R": len(wR), "red_R": len(freely_reduce(wR)),
                "len_1": len(w1), "red_1": len(freely_reduce(w1))}
        for k, v in vals.items():
            if k not in best or v < best[k][0]:
                best[k] = (v, c)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=30)
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)
    out = {"meta": {"experiment": EXPERIMENT_VERSION, "construction": su.CONSTRUCTION_VERSION,
                    "core": CORE_VERSION, "command": " ".join(sys.argv),
                    "python": sys.version.split()[0],
                    "coverage": "exhaustive over the affine family pi(i)=a*i+b, gcd(a,n)=1, "
                                "all n shifts, all words built and executed"},
           "results": {}}
    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        Bn = known.target_diameter(n)
        rec = {"B_n": Bn, "count": 0}
        worst = {}
        for a in range(1, n):
            if gcd(a, n) != 1:
                continue
            for b in range(n):
                pi = tuple((a * i + b) % n for i in range(n))
                rec["count"] += 1
                best = analyze(pi)
                for k, (v, c) in best.items():
                    if k not in worst or v > worst[k][0]:
                        worst[k] = (v, {"a": a, "b": b, "c": c})
        rec["max_over_family_of_min_over_c"] = {
            k: {"value": v, "excess_over_B_n": v - Bn, "witness": w} for k, (v, w) in worst.items()}
        rec["seconds"] = round(time.time() - t0, 1)
        out["results"][n] = rec
        print(f"n={n}: {rec['count']} affine perms, "
              + ", ".join(f"{k} {v['excess_over_B_n']:+d}" for k, v in
                          sorted(rec["max_over_family_of_min_over_c"].items()))
              + f"  ({rec['seconds']} s)", flush=True)
    path = os.path.join(OUT_DIR, f"affine{args.tag}.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=1)
    print("written", path)


if __name__ == "__main__":
    main()
