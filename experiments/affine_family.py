"""Affine family pi(i) = a*i + b (mod n), gcd(a, n) = 1, a != 1, against hypothesis
H10 (CLAIMS C27): the carrier-route variant of construction N1,

    len_R(pi, c) = 2F_c - S_c + R_c,   red_R = |freely reduced word|,

minimised over the shift c, compared with B_n = n(n-1)/2.  Version affine_family-1.0.
Construction constructions/strict_upper.py (strict_upper-1.0), core oracle-1.0;
every word is built, executed with the reference moves (assert in the module)
and freely reduced; the reduced word is re-applied.  d(pi) from certified
tables when available (n <= 10 in git).

Mode "full" (default for n <= --full-max-n): every (a, b) and every c with words.
Mode "best": every (a, b) with statistics only (F_c, S_c, R_c) for all c, then the
words for the worst (a, b) of that n.  Counterexamples (min_c len_R > B_n + 1 or
min_c red_R > B_n) are saved with their words to data/counterexamples/.

Usage: python3 experiments/affine_family.py --min-n 4 --max-n 60 --full-max-n 40
Output: data/runs/h10_verdict/affine_family.json (+ .log)
"""

import argparse
import json
import math
import os
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
for sub in ("oracle", "exact", "constructions"):
    sys.path.insert(0, os.path.join(ROOT, sub))

from moves import apply_word, freely_reduce, identity, CORE_VERSION  # noqa: E402
from bfs import rank_perm, factorials  # noqa: E402
import strict_upper as su  # noqa: E402

EXPERIMENT_VERSION = "affine_family-1.0"
OUT_DIR = os.path.join(ROOT, "data", "runs", "h10_verdict")
CE_DIR = os.path.join(ROOT, "data", "counterexamples")

_TABLES = {}


def dist(p):
    n = len(p)
    if n not in _TABLES:
        path = os.path.join(ROOT, "data", "tables", f"dist_n{n}.bin")
        _TABLES[n] = (open(path, "rb").read(), factorials(n)) if os.path.exists(path) else None
    t = _TABLES[n]
    return None if t is None else t[0][rank_perm(tuple(p), t[1])]


def B(n):
    return n * (n - 1) // 2


def affine(a, b, n):
    return tuple((a * i + b) % n for i in range(n))


def stats_len(pi, c, n):
    f = [(pi[i] + c) % n for i in range(n)]
    F = S = 0
    carriers = []
    K = 0
    for cyc in su.cycles_of(f):
        cd = su.cycle_data(cyc, n)
        F += cd["F"]
        S += cd["S"]
        K += 1
        carriers.append(cyc[su.choose_carrier(cyc, cd, n)])
    R = len(su.carrier_route(carriers, c, n)[1])
    return {"c": c, "K": K, "F": F, "S": S, "R": R, "len_R": 2 * F - S + R}


def with_words(pi, n):
    """All shifts with words: per-c table incl. reduced length; words kept."""
    rows = []
    for c in range(n):
        sw = su.shift_word(pi, c, "carrier")
        st = sw["stats"]
        red = freely_reduce(sw["word"])
        assert apply_word(pi, red) == identity(n)
        rows.append({"c": c, "K": st["n_cycles"], "F": st["F_c"], "S": st["S_c"], "R": st["route_len"],
                     "len_R": st["len"], "red_R": len(red), "N_X": st["N_X"], "N_rot": st["N_rot"],
                     "carriers": st["carriers"], "cycle_lengths": sorted(lw["cd"]["k"] for lw in sw["locals"]),
                     "word": sw["word"], "red_word": red})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=60)
    ap.add_argument("--full-max-n", type=int, default=40)
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(CE_DIR, exist_ok=True)
    t0 = time.time()
    results = {}
    counterexamples = []
    for n in range(args.min_n, args.max_n + 1):
        tn = time.time()
        full = n <= args.full_max_n
        rec = {"n": n, "B_n": B(n), "mode": "full" if full else "best", "pairs": 0,
               "worst_len": None, "worst_red": None, "hist_min_len_R_minus_B": {}, "hist_min_red_R_minus_B": {},
               "worst_red_minus_d": None, "min_B_minus_d": None}
        hist_len, hist_red = {}, {}
        for a in range(2, n):
            if math.gcd(a, n) != 1:
                continue
            for b in range(n):
                pi = affine(a, b, n)
                rec["pairs"] += 1
                if full:
                    rows = with_words(pi, n)
                    ml = min(r["len_R"] for r in rows) - B(n)
                    mr = min(r["red_R"] for r in rows) - B(n)
                    hist_red[mr] = hist_red.get(mr, 0) + 1
                    if rec["worst_red"] is None or mr > rec["worst_red"]["min_red_R_minus_B"]:
                        rec["worst_red"] = {"a": a, "b": b, "min_red_R_minus_B": mr, "min_len_R_minus_B": ml}
                    d = dist(pi)
                    if d is not None:
                        loss = min(r["red_R"] for r in rows) - d
                        if rec["worst_red_minus_d"] is None or loss > rec["worst_red_minus_d"]["loss"]:
                            rec["worst_red_minus_d"] = {"a": a, "b": b, "d": d, "min_red_R": min(r["red_R"] for r in rows), "loss": loss}
                        if rec["min_B_minus_d"] is None or B(n) - d < rec["min_B_minus_d"]["B_minus_d"]:
                            rec["min_B_minus_d"] = {"a": a, "b": b, "d": d, "B_minus_d": B(n) - d}
                    if ml > 1 or mr > 0:
                        counterexamples.append({"n": n, "a": a, "b": b, "pi": list(pi), "min_len_R_minus_B": ml,
                                                "min_red_R_minus_B": mr, "d": dist(pi), "per_c": rows})
                else:
                    ml = min(stats_len(pi, c, n)["len_R"] for c in range(n)) - B(n)
                hist_len[ml] = hist_len.get(ml, 0) + 1
                if rec["worst_len"] is None or ml > rec["worst_len"]["min_len_R_minus_B"]:
                    rec["worst_len"] = {"a": a, "b": b, "min_len_R_minus_B": ml}
        if not full:
            a, b = rec["worst_len"]["a"], rec["worst_len"]["b"]
            rows = with_words(affine(a, b, n), n)
            ml = min(r["len_R"] for r in rows) - B(n)
            mr = min(r["red_R"] for r in rows) - B(n)
            assert ml == rec["worst_len"]["min_len_R_minus_B"]
            rec["worst_len"]["min_red_R_minus_B"] = mr
            if ml > 1 or mr > 0:
                counterexamples.append({"n": n, "a": a, "b": b, "pi": list(affine(a, b, n)), "min_len_R_minus_B": ml,
                                        "min_red_R_minus_B": mr, "d": None, "per_c": rows})
        rec["hist_min_len_R_minus_B"] = {str(k): v for k, v in sorted(hist_len.items())}
        rec["hist_min_red_R_minus_B"] = {str(k): v for k, v in sorted(hist_red.items())}
        rec["seconds"] = round(time.time() - tn, 1)
        results[n] = rec
        print(f"n={n} ({rec['mode']}, {rec['pairs']} pairs, {rec['seconds']}s): worst min_c len_R - B = {rec['worst_len']}; "
              f"worst red = {rec['worst_red']}; worst red - d = {rec['worst_red_minus_d']}; min B - d = {rec['min_B_minus_d']}", flush=True)
    # save counterexamples: the first for each n (words), others without words
    for ce in counterexamples:
        path = os.path.join(CE_DIR, f"H10_affine_n{ce['n']}_a{ce['a']}_b{ce['b']}.json")
        json.dump({"meta": {"experiment": EXPERIMENT_VERSION, "construction": su.CONSTRUCTION_VERSION, "core": CORE_VERSION,
                            "claim": "H10 / C27: min_c [2F_c - S_c + R_c] <= B_n + 1 and min_c red_R <= B_n",
                            "input": f"pi(i) = {ce['a']}*i + {ce['b']} mod {ce['n']}"}, **ce}, open(path, "w"), indent=1)
    report = {"meta": {"experiment": EXPERIMENT_VERSION, "construction": su.CONSTRUCTION_VERSION, "core": CORE_VERSION,
                       "command": " ".join(sys.argv), "python": sys.version.split()[0], "seconds": round(time.time() - t0, 1)},
              "coverage": {"full": f"{args.min_n}<=n<={min(args.max_n, args.full_max_n)}: all (a,b), all c, words built/executed/reduced",
                           "best": f"{args.full_max_n + 1}<=n<={args.max_n}: all (a,b) statistics for all c; words for the worst (a,b)",
                           "d": "certified tables n<=10"},
              "results": results,
              "counterexamples_saved": [os.path.basename(os.path.join(CE_DIR, f"H10_affine_n{ce['n']}_a{ce['a']}_b{ce['b']}.json")) for ce in counterexamples]}
    json.dump(report, open(os.path.join(OUT_DIR, "affine_family.json"), "w"), indent=1)
    print("DONE", report["meta"]["seconds"], "s;", len(counterexamples), "counterexamples saved")


if __name__ == "__main__":
    main()
