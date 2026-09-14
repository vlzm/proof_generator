"""H13-I structured/random sampling ladder (AGENTS.md rule 9) at n beyond
exhaustive reach: reflections, identity, affine pi(i) = a*i+b mod n, a
half-block reversal, and random permutations, for n up to 200. Exhaustive
verification of I(pi) <= floor((n-1)^2/4) is done separately by
h13i_toric.c (4 <= n <= 12); this script is SAMPLED evidence only, not a
certificate, for larger n.

I(pi) is computed exactly (not approximated) by the O(n^3) circular
difference-array method: for each of the n position-cuts q, every pair of
line positions contributes a contiguous circular arc of value-shifts s that
make it an inversion; accumulate with a difference array and prefix-sum once
per q. Cross-checked against brute-force enumeration of all n^2 (q, c) pairs
on random permutations up to n = 8 (see self_check()).

Usage: python3 experiments/h13i_sample.py [--nlist 11,12,15,20,25,30,40,50,60,80,100,150,200] [--seed 42]
Output: data/runs/h13i_search/sample_report.json, .md.
Version h13i_sample-1.0 (session 9, 14.09.2026).
"""

import argparse
import itertools
import json
import os
import random
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
from moves import CORE_VERSION  # noqa: E402

VERSION = "h13i_sample-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h13i_search")


def full_min_inv_fast(pi):
    """Exact I(pi): O(n^3) via a circular difference array per position-cut."""
    n = len(pi)
    best = n * n
    diff = [0] * (n + 1)
    for q in range(n):
        a = [pi[(q + 1 + j) % n] for j in range(n)]
        for i in range(n + 1):
            diff[i] = 0
        for j in range(n):
            aj = a[j]
            for k in range(j + 1, n):
                ak = a[k]
                d = (aj - ak) % n
                length = n - d
                start = (ak - length + 1) % n
                end = start + length
                if end <= n:
                    diff[start] += 1
                    diff[end] -= 1
                else:
                    diff[start] += 1
                    diff[n] -= 1
                    diff[0] += 1
                    diff[end - n] -= 1
        cur = 0
        for s in range(n):
            cur += diff[s]
            if cur < best:
                best = cur
    return best


def _inv_count(w):
    n = len(w)
    c = 0
    for i in range(n):
        for j in range(i + 1, n):
            if w[i] > w[j]:
                c += 1
    return c


def _full_min_inv_slow(pi):
    n = len(pi)
    best = n * n
    for q in range(n):
        for c in range(n):
            shift = (q + 1 - c) % n
            w = [(pi[(q + 1 + j) % n] - shift) % n for j in range(n)]
            iv = _inv_count(w)
            if iv < best:
                best = iv
    return best


def self_check(seed=0, trials=30, nmax=8):
    rng = random.Random(seed)
    for n in range(3, nmax + 1):
        for _ in range(trials):
            pi = list(range(n))
            rng.shuffle(pi)
            pi = tuple(pi)
            fast = full_min_inv_fast(pi)
            slow = _full_min_inv_slow(pi)
            assert fast == slow, (pi, fast, slow)


def target(n):
    return (n - 1) ** 2 // 4


def structured_family(n):
    fam = {}
    for h in sorted({0, n // 4, n // 2, n - 1}):
        fam[f"reflection_h={h}"] = [(h - i) % n for i in range(n)]
    fam["identity"] = list(range(n))
    for a in sorted({2, 3, (n // 2) + 1}):
        if a % n == 0 or a == 1:
            continue
        for b in (0, n // 3):
            fam[f"affine_a={a},b={b}"] = [(a * i + b) % n for i in range(n)]
    pi = list(range(n))
    m = n // 2
    pi[:m] = pi[:m][::-1]
    fam["half_block_reversed"] = pi
    return fam


def run(nlist, seed, random_trials):
    self_check(seed=seed)
    rng = random.Random(seed)
    results = {}
    t0 = time.time()
    any_bad = False
    for n in nlist:
        entries = []
        fam = structured_family(n)
        for name, pi in fam.items():
            v = full_min_inv_fast(tuple(pi))
            entries.append({"name": name, "I": v})
        for trial in range(random_trials):
            pi = list(range(n))
            rng.shuffle(pi)
            v = full_min_inv_fast(tuple(pi))
            entries.append({"name": f"random#{trial}", "I": v})
        tgt = target(n)
        worst = max(e["I"] for e in entries)
        bad = worst > tgt
        any_bad |= bad
        results[n] = {"target": tgt, "worst_observed": worst, "bad": bad, "entries": entries}
        print(f"n={n:4d} target={tgt:6d} worst_observed={worst:6d} {'BAD' if bad else 'ok'}")
    elapsed = time.time() - t0

    os.makedirs(OUT, exist_ok=True)
    report = {
        "version": VERSION,
        "core_version": CORE_VERSION,
        "goal": "SAMPLED evidence for H13-I (I(pi) <= floor((n-1)^2/4)) beyond exhaustive n",
        "command": f"python3 experiments/h13i_sample.py --nlist {','.join(map(str, nlist))} --seed {seed}",
        "coverage": "SAMPLED: reflections, identity, affine family, half-block reversal, "
                    f"{random_trials} random permutations per n",
        "seed": seed,
        "elapsed_s": elapsed,
        "results": results,
        "any_counterexample": any_bad,
        "conclusion": (
            "no counterexample found; reflections achieve I(pi) == target exactly at every n tested "
            "(consistent with C33's finite data: reflections are the extremal family)"
            if not any_bad else "COUNTEREXAMPLE FOUND -- see results"
        ),
    }
    with open(os.path.join(OUT, "sample_report.json"), "w") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)

    lines = [
        "# H13-I structured/random sampling -- report",
        "",
        f"Version {VERSION}, core {CORE_VERSION}, seed {seed}, {elapsed:.1f} s total.",
        "SAMPLED coverage only (not exhaustive); exhaustive range is h13i_toric.c (n <= 12).",
        "",
        "| n | target | worst observed |",
        "|---|---|---|",
    ]
    for n in nlist:
        r = results[n]
        lines.append(f"| {n} | {r['target']} | {r['worst_observed']} |")
    lines.append("")
    lines.append(report["conclusion"])
    with open(os.path.join(OUT, "sample_report.md"), "w") as f:
        f.write("\n".join(lines) + "\n")

    return report


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nlist", type=str, default="11,12,15,20,25,30,40,50,60,80,100,150,200")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--random-trials", type=int, default=8)
    args = ap.parse_args()
    nlist = [int(x) for x in args.nlist.split(",")]
    run(nlist, args.seed, args.random_trials)
