"""Evaluate the candidate sweep construction (constructions/candidate_sweep.py)
against the certified tables and against B_n.  Version sweep_eval-1.0.

Usage: python3 experiments/sweep_eval.py --perms 8 [--families 60] [--samples 200] [--workers 4]
Portfolio = min over the N1 carrier-route variant (strict_upper, min over c) and the sweep candidate.
Output: data/runs/sweep_eval/report.json, report.md (+ log)
"""

import argparse
import itertools
import json
import os
import random
import sys
import time
from collections import Counter
from math import gcd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
sys.path.insert(0, os.path.join(ROOT, "exact"))
sys.path.insert(0, os.path.join(ROOT, "constructions"))

from moves import sigma, rev, CORE_VERSION  # noqa: E402
from bfs import factorials, rank_perm  # noqa: E402
import candidate_sweep as cs  # noqa: E402
import strict_upper as su  # noqa: E402

VERSION = "sweep_eval-1.0"
OUT_DIR = os.path.join(ROOT, "data", "runs", "sweep_eval")


def table(n):
    path = os.path.join(ROOT, "data", "tables", f"dist_n{n}.bin")
    if not os.path.exists(path):
        return None
    with open(path, "rb") as f:
        return f.read()


def n1_best(pi):
    n = len(pi)
    return min(su.shift_word(pi, c, route="carrier")["stats"]["len"] for c in range(n))


def _chunk(args):
    n, first = args
    B = n * (n - 1) // 2
    tab = table(n)
    fact = factorials(n)
    exc, red_exc, loss, port, wins = Counter(), Counter(), Counter(), Counter(), Counter()
    worst = []
    cnt = 0
    for rest in itertools.permutations([i for i in range(n) if i != first]):
        pi = (first,) + rest
        cnt += 1
        r = cs.best_sweep(pi)
        e = r["len"] - B
        exc[e] += 1
        red_exc[r["reduced_len"] - B] += 1
        n1 = n1_best(pi)
        pe = min(e, n1 - B)
        port[pe] += 1
        wins["sweep" if r["len"] < n1 else ("tie" if r["len"] == n1 else "n1")] += 1
        d = int(tab[rank_perm(pi, fact)]) if tab is not None else None
        if d is not None:
            loss[r["len"] - d] += 1
        if e >= 0 or pe >= 0:
            worst.append({"pi": list(pi), "len": r["len"], "reduced": r["reduced_len"],
                          "c": r["c"], "dir": r["dir"], "thr": r["thr"], "mode": r["mode"],
                          "passes": r["passes"], "N_X": r["N_X"], "N_rot": r["N_rot"],
                          "word": r["word"], "d": d, "N1_min_len": n1})
    return cnt, exc, red_exc, loss, port, wins, worst


def eval_perms(n, log, workers):
    t0 = time.time()
    if workers > 1:
        from multiprocessing import Pool
        with Pool(workers) as pool:
            parts = pool.map(_chunk, [(n, f) for f in range(n)])
    else:
        parts = [_chunk((n, f)) for f in range(n)]
    cnt = sum(p[0] for p in parts)
    exc, red_exc, loss, port, wins = Counter(), Counter(), Counter(), Counter(), Counter()
    worst = []
    for p in parts:
        exc.update(p[1]); red_exc.update(p[2]); loss.update(p[3]); port.update(p[4]); wins.update(p[5])
        worst.extend(p[6])
    sec = time.time() - t0
    log(f"n={n}: {cnt} perms, {sec:.0f} s; sweep: max len-B_n = {max(exc)} ({exc[max(exc)]} inputs), "
        f"max reduced-B_n = {max(red_exc)}, max len-d = {max(loss) if loss else None}; "
        f"portfolio min(N1_R, sweep): max -B_n = {max(port)} ({port[max(port)]} inputs); wins {dict(wins)}")
    return {"n": n, "count": cnt, "seconds": round(sec, 1),
            "hist_len_minus_B": dict(sorted(exc.items())),
            "hist_reduced_minus_B": dict(sorted(red_exc.items())),
            "hist_len_minus_d": dict(sorted(loss.items())),
            "hist_portfolio_minus_B": dict(sorted(port.items())),
            "wins": dict(wins),
            "worst": sorted(worst, key=lambda w: (-min(w["len"], w["N1_min_len"]), -w["len"]))[:60]}


def eval_families(max_n, log, samples, seed):
    rng = random.Random(seed)
    rows = []
    t0 = time.time()
    for n in range(4, max_n + 1):
        B = n * (n - 1) // 2
        fam = {}
        refl = [tuple((h - i) % n for i in range(n)) for h in range(n)]
        fam["reflections"] = max(cs.best_sweep(p)["len"] - B for p in refl)
        fam["sigma_n"] = cs.best_sweep(sigma(n))["len"] - B
        fam["rev_n"] = cs.best_sweep(rev(n))["len"] - B
        ma, arg, mp, argp = None, None, None, None
        for a in range(2, n - 1):
            if gcd(a, n) != 1:
                continue
            for b in (range(n) if n <= 20 else (0, 1, n // 2, n - 1)):
                pi = tuple((a * i + b) % n for i in range(n))
                v = cs.best_sweep(pi)["len"] - B
                if ma is None or v > ma:
                    ma, arg = v, (a, b)
                pv = min(v, n1_best(pi) - B)
                if mp is None or pv > mp:
                    mp, argp = pv, (a, b)
        fam["affine_no_refl"] = ma
        fam["affine_b_coverage"] = "all b" if n <= 20 else "b in {0, 1, n//2, n-1}"
        fam["affine_argmax"] = arg
        fam["affine_portfolio"] = mp
        fam["affine_portfolio_argmax"] = argp
        fam["rotations_sigma"] = max(cs.best_sweep(tuple(sigma(n)[(i + k) % n] for i in range(n)))["len"] - B
                                     for k in range(n))
        if samples:
            mx, mxp = None, None
            for _ in range(samples):
                p = list(range(n))
                rng.shuffle(p)
                v = cs.best_sweep(tuple(p))["len"] - B
                pv = min(v, n1_best(tuple(p)) - B)
                mx = v if mx is None else max(mx, v)
                mxp = pv if mxp is None else max(mxp, pv)
            fam["random"] = mx
            fam["random_portfolio"] = mxp
        rows.append({"n": n, **fam})
        log(f"families n={n}: {fam}")
    log(f"families total {time.time() - t0:.0f} s")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perms", type=int, default=8)
    ap.add_argument("--families", type=int, default=0)
    ap.add_argument("--samples", type=int, default=0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=1)
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)
    logf = open(os.path.join(OUT_DIR, "report.log"), "a")

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + "\n")
        logf.flush()

    log(f"== {VERSION} {cs.CANDIDATE_VERSION} {CORE_VERSION} args={vars(args)}")
    out = {"meta": {"version": VERSION, "candidate": cs.CANDIDATE_VERSION, "core": CORE_VERSION,
                    "construction_ref": su.CONSTRUCTION_VERSION, "args": vars(args)},
           "perms": [], "families": []}
    for n in range(4, args.perms + 1):
        out["perms"].append(eval_perms(n, log, args.workers))
    if args.families:
        out["families"] = eval_families(args.families, log, args.samples, args.seed)
    path = os.path.join(OUT_DIR, f"report_perms{args.perms}_fam{args.families}.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=1)
    log(f"written {path}")


if __name__ == "__main__":
    main()
