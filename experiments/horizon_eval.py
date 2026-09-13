"""Evaluate the "horizon" head policy (candidate_sweep sweep-1.2, mode "horizon")
against the certified tables, B_n and the N1 carrier-route variant.  H12 step (1).

Usage:
  python3 experiments/horizon_eval.py --perms 9 --workers 4 --horizon 3 --weight 1
  python3 experiments/horizon_eval.py --families 20 --samples 30 --horizon 3
  python3 experiments/horizon_eval.py --large 24,30,40 --samples 10 --workers 4 --horizon 3
Modes evaluated per input: "horizon" alone (the candidate policy), and the
portfolio min(N1_R over c, horizon).  For --perms also len - d from the table.
Output: data/runs/horizon_eval/report_<tag>.json, report.log
Version horizon_eval-1.0.
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

VERSION = "horizon_eval-1.1"
OUT_DIR = os.path.join(ROOT, "data", "runs", "horizon_eval")
CFG = {"horizon": 3, "weight": 1.0, "thrs": (1, 2), "fallback": False}


def table(n):
    path = os.path.join(ROOT, "data", "tables", f"dist_n{n}.bin")
    if not os.path.exists(path):
        return None
    with open(path, "rb") as f:
        return f.read()


def n1_best(pi):
    n = len(pi)
    return min(su.shift_word(pi, c, route="carrier")["stats"]["len"] for c in range(n))


def hz(pi):
    return cs.best_sweep(pi, thrs=CFG["thrs"], modes=("horizon",),
                         horizon=CFG["horizon"], weight=CFG["weight"], fallback=CFG["fallback"])


def _chunk(args):
    n, first, cfg = args
    CFG.update(cfg)
    B = n * (n - 1) // 2
    tab = table(n)
    fact = factorials(n)
    exc, loss, port, wins = Counter(), Counter(), Counter(), Counter()
    worst = []
    cnt = 0
    for rest in itertools.permutations([i for i in range(n) if i != first]):
        pi = (first,) + rest
        cnt += 1
        r = hz(pi)
        e = r["len"] - B
        exc[e] += 1
        d = int(tab[rank_perm(pi, fact)]) if tab is not None else None
        if d is not None:
            loss[r["len"] - d] += 1
        n1 = n1_best(pi) if e >= -1 else None
        if n1 is not None:
            pe = min(e, n1 - B)
            port[pe] += 1
            wins["horizon" if r["len"] < n1 else ("tie" if r["len"] == n1 else "n1")] += 1
        if e >= 0:
            worst.append({"pi": list(pi), "len": r["len"], "c": r["c"], "dir": r["dir"],
                          "thr": r["thr"], "N_X": r["N_X"], "N_rot": r["N_rot"],
                          "word": r["word"], "d": d, "N1_min_len": n1})
    return cnt, exc, loss, port, wins, worst


def eval_perms(n, log, workers):
    t0 = time.time()
    jobs = [(n, f, dict(CFG)) for f in range(n)]
    if workers > 1:
        from multiprocessing import Pool
        with Pool(workers) as pool:
            parts = pool.map(_chunk, jobs)
    else:
        parts = [_chunk(j) for j in jobs]
    cnt = sum(p[0] for p in parts)
    exc, loss, port, wins = Counter(), Counter(), Counter(), Counter()
    worst = []
    for p in parts:
        exc.update(p[1]); loss.update(p[2]); port.update(p[3]); wins.update(p[4])
        worst.extend(p[5])
    sec = time.time() - t0
    log(f"n={n}: {cnt} perms, {sec:.0f} s; horizon: max len-B_n = {max(exc)} ({exc[max(exc)]} inputs), "
        f"max len-d = {max(loss) if loss else None}; portfolio (only inputs with len >= B_n-1): "
        f"max -B_n = {max(port) if port else None}; wins there {dict(wins)}")
    return {"n": n, "count": cnt, "seconds": round(sec, 1),
            "hist_len_minus_B": dict(sorted(exc.items())),
            "hist_len_minus_d": dict(sorted(loss.items())),
            "hist_portfolio_minus_B_on_len_ge_Bminus1": dict(sorted(port.items())),
            "wins_on_len_ge_Bminus1": dict(wins),
            "worst": sorted(worst, key=lambda w: -w["len"])[:60]}


def eval_families(max_n, log, samples, seed):
    rng = random.Random(seed)
    rows = []
    t0 = time.time()
    for n in range(4, max_n + 1):
        B = n * (n - 1) // 2
        fam = {}
        refl = [tuple((h - i) % n for i in range(n)) for h in range(n)]
        rv = [hz(p)["len"] - B for p in refl]
        fam["reflections"] = max(rv)
        fam["reflections_argmax_h"] = rv.index(max(rv))
        fam["sigma_n"] = hz(sigma(n))["len"] - B
        fam["rev_n"] = hz(rev(n))["len"] - B
        ma, arg = None, None
        for a in range(2, n - 1):
            if gcd(a, n) != 1:
                continue
            for b in (range(n) if n <= 20 else (0, 1, n // 2, n - 1)):
                pi = tuple((a * i + b) % n for i in range(n))
                v = hz(pi)["len"] - B
                if ma is None or v > ma:
                    ma, arg = v, (a, b)
        fam["affine_no_refl"] = ma
        fam["affine_b_coverage"] = "all b" if n <= 20 else "b in {0, 1, n//2, n-1}"
        fam["affine_argmax"] = arg
        fam["rotations_sigma"] = max(hz(tuple(sigma(n)[(i + k) % n] for i in range(n)))["len"] - B
                                     for k in range(n))
        # transposition products / short cycles
        fam["adjacent_transpositions"] = max(
            hz(tuple((i + 1 if i % 2 == 0 and i + 1 < n else (i - 1 if i % 2 == 1 else i)) for i in range(n)))["len"] - B,
            hz(tuple(range(n)))["len"] - B)
        if samples:
            mx = None
            for _ in range(samples):
                p = list(range(n))
                rng.shuffle(p)
                v = hz(tuple(p))["len"] - B
                mx = v if mx is None else max(mx, v)
            fam["random"] = mx
        rows.append({"n": n, **fam})
        log(f"families n={n}: {fam} ({time.time() - t0:.0f} s)")
    return rows


def _one(args):
    name, pi, cfg = args
    CFG.update(cfg)
    n = len(pi)
    B = n * (n - 1) // 2
    r = hz(pi)
    return name, r["len"] - B, r["N_X"], r["N_rot"], r["c"]


def eval_large(ns, log, samples, seed, workers):
    """Families at selected large n, parallel over inputs; full min over c."""
    from multiprocessing import Pool
    rng = random.Random(seed)
    rows = []
    for n in ns:
        t0 = time.time()
        inputs = []
        for h in range(n):
            inputs.append((f"refl_h{h}", tuple((h - i) % n for i in range(n))))
        inputs.append(("rev_n", rev(n)))
        aa = [a for a in range(2, n - 1) if gcd(a, n) == 1]
        if n > 40:
            aa = aa[:3] + aa[-3:] + [a for a in aa if abs(a - n // 2) <= 2]
        for a in sorted(set(aa)):
            for b in (0, 1, n // 2, n - 1):
                inputs.append((f"affine_{a}_{b}", tuple((a * i + b) % n for i in range(n))))
        for k in range(samples):
            p = list(range(n))
            rng.shuffle(p)
            inputs.append((f"random_{k}", tuple(p)))
        jobs = [(name, pi, dict(CFG)) for name, pi in inputs]
        with Pool(workers) as pool:
            res = pool.map(_one, jobs, chunksize=1)
        fam = {}
        for name, e, nx, nr, c in res:
            key = name.split("_")[0]
            if key not in fam or e > fam[key][0]:
                fam[key] = (e, name, nx, nr, c)
        row = {"n": n, "inputs": len(inputs), "seconds": round(time.time() - t0),
               "max_excess": {k: v[0] for k, v in fam.items()},
               "argmax": {k: v[1] for k, v in fam.items()},
               "all": [{"name": r[0], "excess": r[1], "N_X": r[2], "N_rot": r[3], "c": r[4]} for r in res]}
        rows.append(row)
        log(f"large n={n}: {len(inputs)} inputs, {row['seconds']} s; max excess {row['max_excess']} argmax {row['argmax']}")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perms", type=int, default=0)
    ap.add_argument("--perms_from", type=int, default=4)
    ap.add_argument("--families", type=int, default=0)
    ap.add_argument("--samples", type=int, default=0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--horizon", type=int, default=3)
    ap.add_argument("--weight", type=float, default=1.0)
    ap.add_argument("--thrs", type=str, default="1,2")
    ap.add_argument("--tag", type=str, default="")
    ap.add_argument("--fallback", action="store_true", help="long-arc fallback on deadlock (sweep-1.3)")
    ap.add_argument("--large", type=str, default="", help="comma list of n for family checks in parallel")
    args = ap.parse_args()
    CFG["horizon"] = args.horizon
    CFG["weight"] = args.weight
    CFG["thrs"] = tuple(int(t) for t in args.thrs.split(","))
    CFG["fallback"] = bool(args.fallback)
    os.makedirs(OUT_DIR, exist_ok=True)
    logf = open(os.path.join(OUT_DIR, "report.log"), "a")

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + "\n")
        logf.flush()

    log(f"== {VERSION} {cs.CANDIDATE_VERSION} {CORE_VERSION} args={vars(args)}")
    out = {"meta": {"version": VERSION, "candidate": cs.CANDIDATE_VERSION, "core": CORE_VERSION,
                    "construction_ref": su.CONSTRUCTION_VERSION, "args": vars(args), "cfg": dict(CFG)},
           "perms": [], "families": [], "large": []}
    for n in range(args.perms_from, args.perms + 1):
        out["perms"].append(eval_perms(n, log, args.workers))
    if args.families:
        out["families"] = eval_families(args.families, log, args.samples, args.seed)
    if args.large:
        out["large"] = eval_large([int(x) for x in args.large.split(",")], log, args.samples, args.seed, args.workers)
    tag = args.tag or f"H{args.horizon}w{args.weight}_perms{args.perms}_fam{args.families}"
    path = os.path.join(OUT_DIR, f"report_{tag}.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=1)
    log(f"written {path}")


if __name__ == "__main__":
    main()
