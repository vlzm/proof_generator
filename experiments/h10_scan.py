"""Exhaustive scan of hypothesis H10 (CLAIMS C27) for one n: for every pi in S_n
and every shift c the carrier-route variant of construction N1 is evaluated,

    len_R(pi, c) = 2F_c - S_c + R_c,      red_R = |freely reduced word|,

and the minima over c are compared with B_n = n(n-1)/2.  Version h10_scan-1.0.
Construction: constructions/strict_upper.py (strict_upper-1.0), core oracle-1.0.

Per permutation: statistics (F_c, S_c, carriers, R_c) for ALL n shifts without
building words; the full word is built, executed (assert: sorts pi) and freely
reduced for every shift with len_R >= B_n - SLACK (default SLACK = 0), i.e.
exactly where the reduced length matters for the second part of H10, and always
for the best shift.  Records: histograms of min_c len_R - B_n and
min_c red_R - B_n, every pi with min_c len_R >= B_n (with its per-c table),
counts of words built, time; checkpoint every CHUNK permutations (resumable).

Usage: python3 experiments/h10_scan.py --n 10 [--workers 4] [--slack 0]
Output: data/runs/h10_scan/scan_n<n>.json (+ .log), checkpoint scan_n<n>.ckpt.json
"""

import argparse
import itertools
import json
import os
import sys
import time
import collections
from multiprocessing import Pool

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
for sub in ("oracle", "constructions"):
    sys.path.insert(0, os.path.join(ROOT, sub))

from moves import apply_word, freely_reduce, identity, CORE_VERSION  # noqa: E402
import strict_upper as su  # noqa: E402

SCAN_VERSION = "h10_scan-1.0"
OUT_DIR = os.path.join(ROOT, "data", "runs", "h10_scan")


def B(n):
    return n * (n - 1) // 2


def shift_stats(pi, c, n):
    f = [(pi[i] + c) % n for i in range(n)]
    F = S = 0
    carriers = []
    for cyc in su.cycles_of(f):
        cd = su.cycle_data(cyc, n)
        F += cd["F"]
        S += cd["S"]
        carriers.append(cyc[su.choose_carrier(cyc, cd, n)])
    pos, letters = su.carrier_route(carriers, c, n)
    return 2 * F - S + len(letters), F, S, len(letters), carriers


def scan_perm(pi, n, slack):
    """Returns (min_c len_R - B_n, min_c red_R - B_n, words_built, table or None)."""
    Bn = B(n)
    stats = [shift_stats(pi, c, n) for c in range(n)]
    best_len = min(s[0] for s in stats)
    best_c = min(range(n), key=lambda c: (stats[c][0], c))
    best_red = None
    built = 0
    for c in range(n):
        if stats[c][0] >= Bn - slack or c == best_c:
            sw = su.shift_word(pi, c, "carrier")          # asserts the word sorts pi
            assert len(sw["word"]) == stats[c][0], (pi, c)
            red = len(freely_reduce(sw["word"]))
            built += 1
            if best_red is None or red < best_red:
                best_red = red
    table = None
    if best_len >= Bn:
        table = [{"c": c, "len_R": s[0], "F": s[1], "S": s[2], "R": s[3], "carriers": s[4]} for c, s in enumerate(stats)]
    return best_len - Bn, best_red - Bn, built, table


def worker(args):
    n, first, slack = args
    hist_len = collections.Counter()
    hist_red = collections.Counter()
    tight = []
    built = 0
    count = 0
    rest = [i for i in range(n) if i != first]
    for tail in itertools.permutations(rest):
        pi = (first,) + tail
        a, b, w, table = scan_perm(pi, n, slack)
        hist_len[a] += 1
        hist_red[b] += 1
        built += w
        count += 1
        if table is not None:
            tight.append({"pi": list(pi), "min_len_R_minus_B": a, "min_red_R_minus_B": b, "table": table})
    return first, count, dict(hist_len), dict(hist_red), tight, built


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--slack", type=int, default=0)
    args = ap.parse_args()
    n = args.n
    os.makedirs(OUT_DIR, exist_ok=True)
    ckpt_path = os.path.join(OUT_DIR, f"scan_n{n}.ckpt.json")
    out_path = os.path.join(OUT_DIR, f"scan_n{n}.json")
    done = {}
    if os.path.exists(ckpt_path):
        done = {int(k): v for k, v in json.load(open(ckpt_path))["chunks"].items()}
        print(f"resuming: chunks done {sorted(done)}", flush=True)
    todo = [(n, first, args.slack) for first in range(n) if first not in done]
    t0 = time.time()
    with Pool(args.workers) as pool:
        for first, count, hl, hr, tight, built in pool.imap_unordered(worker, todo):
            done[first] = {"count": count, "hist_len": hl, "hist_red": hr, "tight": tight, "built": built,
                           "elapsed": time.time() - t0}
            json.dump({"chunks": {str(k): v for k, v in done.items()}}, open(ckpt_path, "w"))
            print(f"chunk pi[0]={first}: {count} perms, hist len_R-B {sorted((int(k), v) for k, v in hl.items())[-4:]}, "
                  f"hist red_R-B {sorted((int(k), v) for k, v in hr.items())[-3:]}, tight {len(tight)}, "
                  f"words {built}, {time.time() - t0:.0f}s", flush=True)
    hist_len = collections.Counter()
    hist_red = collections.Counter()
    tight = []
    total = built = 0
    for v in done.values():
        hist_len.update({int(k): c for k, c in v["hist_len"].items()})
        hist_red.update({int(k): c for k, c in v["hist_red"].items()})
        tight.extend(v["tight"])
        total += v["count"]
        built += v["built"]
    report = {
        "meta": {"scan": SCAN_VERSION, "construction": su.CONSTRUCTION_VERSION, "core": CORE_VERSION,
                 "command": " ".join(sys.argv), "n": n, "B_n": B(n), "slack": args.slack,
                 "workers": args.workers, "python": sys.version.split()[0]},
        "coverage": {"permutations": total, "shifts_stats": total * n, "words_built_executed": built,
                     "note": "stats for all shifts; words built, executed and freely reduced for the best shift "
                             "and for every shift with len_R >= B_n - slack"},
        "hist_min_len_R_minus_B_n": {str(k): v for k, v in sorted(hist_len.items())},
        "hist_min_red_R_minus_B_n": {str(k): v for k, v in sorted(hist_red.items())},
        "max_min_len_R_minus_B_n": max(hist_len),
        "max_min_red_R_minus_B_n": max(hist_red),
        "tight_inputs": tight,
        "seconds": time.time() - t0,
    }
    json.dump(report, open(out_path, "w"), indent=1)
    print(f"DONE n={n}: {total} perms; max min_c len_R - B_n = {max(hist_len)}; "
          f"max min_c red_R - B_n = {max(hist_red)}; tight {len(tight)}; {report['seconds']:.0f}s -> {out_path}", flush=True)


if __name__ == "__main__":
    main()
