"""Evaluate the candidate "frontier" construction (constructions/candidate_frontier.py,
mode "weighted") against the certified tables and against B_n -- the H12 step 1
criterion (PLAN.md row H12, docs/notes/geodesic_structure.md section 5): a single
head policy with no portfolio, <= B_n on all pi for 4 <= n <= 9 exhaustively and
on reflections/sigma_n/rev_n/affine families for n <= 60.

Usage: python3 experiments/frontier_eval.py --perms 9 [--families 60] [--samples 200]
       [--workers 4]
Also reports the N1 carrier-route variant and the sweep-1.1 candidate side by side,
and the "farthest" mode of candidate_frontier (empirically == sweep's "nearest",
kept only for the negative-result record, see candidate_frontier.py docstring).
Output: data/runs/frontier_eval/report_perms<N>_fam<M>.json, report.md (+ log)
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
import candidate_frontier as cf  # noqa: E402
import candidate_sweep as cs  # noqa: E402
import strict_upper as su  # noqa: E402

VERSION = "frontier_eval-1.0"
OUT_DIR = os.path.join(ROOT, "data", "runs", "frontier_eval")


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
    exc, loss, port = Counter(), Counter(), Counter()
    wins = Counter()
    worst = []
    cnt = 0
    for rest in itertools.permutations([i for i in range(n) if i != first]):
        pi = (first,) + rest
        cnt += 1
        r = cf.best_weighted(pi)
        e = r["len"] - B
        exc[e] += 1
        n1 = n1_best(pi)
        sw = cs.best_sweep(pi)["len"]
        pe = min(e, n1 - B, sw - B)
        port[pe] += 1
        wins["frontier" if r["len"] < min(n1, sw) else ("tie" if r["len"] == min(n1, sw) else "other")] += 1
        d = int(tab[rank_perm(pi, fact)]) if tab is not None else None
        if d is not None:
            loss[r["len"] - d] += 1
        if e > 0 or pe > 0:
            worst.append({"pi": list(pi), "len": r["len"], "c": r["c"], "thr": r["thr"],
                          "tie_pref": r["tie_pref"], "init_dir": r["init_dir"],
                          "word": r["word"], "d": d, "N1_min_len": n1, "sweep_len": sw})
    return cnt, exc, loss, port, wins, worst


def eval_perms(n, log, workers):
    t0 = time.time()
    if workers > 1:
        from multiprocessing import Pool
        with Pool(workers) as pool:
            parts = pool.map(_chunk, [(n, f) for f in range(n)])
    else:
        parts = [_chunk((n, f)) for f in range(n)]
    cnt = sum(p[0] for p in parts)
    exc, loss, port, wins = Counter(), Counter(), Counter(), Counter()
    worst = []
    for p in parts:
        exc.update(p[1]); loss.update(p[2]); port.update(p[3]); wins.update(p[4])
        worst.extend(p[5])
    sec = time.time() - t0
    log(f"n={n}: {cnt} perms, {sec:.0f} s; frontier: max len-B_n = {max(exc)} "
        f"({exc[max(exc)]} inputs), max len-d = {max(loss) if loss else None}; "
        f"portfolio min(N1_R, sweep, frontier): max -B_n = {max(port)} ({port[max(port)]} inputs); "
        f"wins {dict(wins)}")
    return {"n": n, "count": cnt, "seconds": round(sec, 1),
            "hist_len_minus_B": dict(sorted(exc.items())),
            "hist_len_minus_d": dict(sorted(loss.items())),
            "hist_portfolio_minus_B": dict(sorted(port.items())),
            "wins": dict(wins),
            "worst": sorted(worst, key=lambda w: (-w["len"] + w["d"] if w["d"] is not None else 0,
                                                   -w["len"]))[:60]}


def _reflection_h_set(n, full):
    """AGENTS rule 9/10: n <= 40 -- all h (matches session 5's sweep_eval
    convention); n > 40 -- an explicit reduced set (the h = 0 case is the one
    found to be hard in this session, so it is always included), stated in
    the report rather than silently dropped."""
    if full:
        return list(range(n))
    hs = {0, 1, 2, 3, n // 4, n // 2 - 1, n // 2, n // 2 + 1, n - 2, n - 1}
    return sorted(h % n for h in hs)


def _affine_a_set(n, full):
    """Same rule for the multiplier a of the affine family; n > 40 samples a
    BOUNDED number (8) of coprime a evenly spread over the range instead of
    all of them (b is already restricted to {0, 1, n//2, n-1} for n > 20, as
    in sweep_eval/session 5) -- keeps the cost per n bounded as n grows."""
    cop = [a for a in range(2, n - 1) if gcd(a, n) == 1]
    if full or len(cop) <= 8:
        return cop
    idx = sorted({round(i * (len(cop) - 1) / 7) for i in range(8)})
    return [cop[i] for i in idx]


def _family_row(args):
    n, samples, seed = args
    rng = random.Random(seed + n)
    B = n * (n - 1) // 2
    full = n <= 40
    fam = {}
    hs = _reflection_h_set(n, full)
    refl_vals = [(h, cf.best_weighted(tuple((h - i) % n for i in range(n)))["len"] - B) for h in hs]
    fam["reflections_max"] = max(v for _, v in refl_vals)
    fam["reflections_argmax_h"] = max(refl_vals, key=lambda hv: hv[1])[0]
    fam["reflections_h_coverage"] = "all h" if full else f"h in {hs}"
    fam["sigma_n"] = cf.best_weighted(sigma(n))["len"] - B
    fam["rev_n"] = cf.best_weighted(rev(n))["len"] - B
    ma, arg = None, None
    a_set = _affine_a_set(n, full)
    b_range = range(n) if n <= 20 else (0, 1, n // 2, n - 1)
    for a in a_set:
        for b in b_range:
            pi = tuple((a * i + b) % n for i in range(n))
            v = cf.best_weighted(pi)["len"] - B
            if ma is None or v > ma:
                ma, arg = v, (a, b)
    fam["affine_no_refl_max"] = ma
    fam["affine_argmax"] = arg
    fam["affine_a_coverage"] = "all coprime a" if full else f"{len(a_set)} coprime a sampled: {a_set}"
    fam["affine_b_coverage"] = "all b" if n <= 20 else "b in {0, 1, n//2, n-1}"
    if full:
        fam["rotations_sigma_max"] = max(
            cf.best_weighted(tuple(sigma(n)[(i + k) % n] for i in range(n)))["len"] - B
            for k in range(n))
        fam["rotations_sigma_coverage"] = "all k"
    else:
        ks = sorted({0, n // 4, n // 2, 3 * n // 4, n - 1})
        fam["rotations_sigma_max"] = max(
            cf.best_weighted(tuple(sigma(n)[(i + k) % n] for i in range(n)))["len"] - B
            for k in ks)
        fam["rotations_sigma_coverage"] = f"k in {ks}"
    if samples:
        mx = None
        for _ in range(samples):
            p = list(range(n))
            rng.shuffle(p)
            v = cf.best_weighted(tuple(p))["len"] - B
            mx = v if mx is None else max(mx, v)
        fam["random_max"] = mx
        fam["random_samples"] = samples
    return n, fam


def eval_families(max_n, log, samples, seed, workers=1):
    """Logs each n as soon as it is ready (imap_unordered) so a crash or a
    kill loses only the not-yet-finished tail, not the whole run; the JSON
    returned is still sorted by n."""
    t0 = time.time()
    ns = list(range(4, max_n + 1))
    pairs = []
    ckpt = os.path.join(OUT_DIR, "families_checkpoint.json")

    def checkpoint():
        tmp = sorted(pairs, key=lambda r: r[0])
        with open(ckpt, "w") as f:
            json.dump([{"n": n, **fam} for n, fam in tmp], f, indent=1)

    if workers > 1:
        from multiprocessing import Pool
        with Pool(workers) as pool:
            for n, fam in pool.imap_unordered(_family_row, [(n, samples, seed) for n in ns], chunksize=1):
                pairs.append((n, fam))
                log(f"families n={n}: {fam}")
                checkpoint()
    else:
        for n in ns:
            n, fam = _family_row((n, samples, seed))
            pairs.append((n, fam))
            log(f"families n={n}: {fam}")
            checkpoint()
    pairs.sort(key=lambda r: r[0])
    rows = [{"n": n, **fam} for n, fam in pairs]
    log(f"families total {time.time() - t0:.0f} s")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perms", type=int, default=9)
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

    log(f"== {VERSION} {cf.CANDIDATE_VERSION} {CORE_VERSION} args={vars(args)}")
    out = {"meta": {"version": VERSION, "candidate": cf.CANDIDATE_VERSION, "core": CORE_VERSION,
                    "n1_construction": su.CONSTRUCTION_VERSION, "sweep_construction": cs.CANDIDATE_VERSION,
                    "args": vars(args)},
           "perms": [], "families": []}
    for n in range(4, args.perms + 1):
        out["perms"].append(eval_perms(n, log, args.workers))
    if args.families:
        out["families"] = eval_families(args.families, log, args.samples, args.seed, args.workers)
    path = os.path.join(OUT_DIR, f"report_perms{args.perms}_fam{args.families}.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=1)
    log(f"written {path}")


if __name__ == "__main__":
    main()
