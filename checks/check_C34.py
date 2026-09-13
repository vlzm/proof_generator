"""Checker for C34 (line model, session 8): two proved lemmas and their tightness.

  Part A (universal line route, docs/notes/h13_line_model.md, Lemma A):
    the shrinking cocktail route (constructions/line_cocktail.py) sorts every
    line w of n positions with N_X = inv(w) and N_rot = B_n - 1 (free final
    head position; the return walk to circle position 0 is appended only for
    the reference check); each word is executed with the reference moves
    (oracle-1.0) on the circle with the cut edge {n-1, 0} never crossed.  All n! lines for 4 <= n <= AMAX.
    Tightness: the fully reversed line needs N_rot >= B_n - 1 (inv = B_n);
    the exhaustive minimum over lines is in data/runs/line_rotation_budget/.
  Part B (antipodal transposition, Lemma B): for even n, tau_n = id_n with
    positions n/2, n/2 + 1 swapped has I(tau_n) = 1 (one inversion at the cut
    q = n - 1, c = 0; not a rotation) and d(tau_n) = n + 1 (certified tables
    for 4 <= n <= 10; the proof covers n >= 6).  Hence every bound of the form
    d(pi) <= 2 I(pi) + g(n) needs g(n) >= n - 1.
  Part C (profile, from data/runs/line_profile/profile_n*.json): for
    4 <= n <= 10, max_pi (d - I) = floor(n^2/4), max_pi I = floor((n-1)^2/4),
    and the values max_pi (d - 2I) are re-read and printed.
Usage: python3 checks/check_C34.py [--amax 8]
Output: data/runs/check_C34/report.json, report.md.
"""

import argparse
import itertools
import json
import os
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
sys.path.insert(0, os.path.join(ROOT, "exact"))
sys.path.insert(0, os.path.join(ROOT, "constructions"))

from moves import apply_word, identity, CORE_VERSION  # noqa: E402
from bfs import factorials, rank_perm  # noqa: E402
import line_cocktail as lc  # noqa: E402

VERSION = "check_C34-1.0"
OUT = os.path.join(ROOT, "data", "runs", "check_C34")


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def part_a(amax, log):
    rows = []
    ok = True
    for n in range(4, amax + 1):
        t0 = time.time()
        B = n * (n - 1) // 2
        cnt = 0
        for w in itertools.permutations(range(n)):
            word, nx, nrot, _, head = lc.cocktail_word(w)
            # sorted circle with the head at line index `head`; then R^head returns the head to 0
            if apply_word(tuple(w), word + "R" * head) != identity(n) or nx != inversions(w) or nrot != B - 1:
                ok = False
                log(f"FAIL part A: n={n} w={w} word={word} N_X={nx} N_rot={nrot}")
                break
            # the head never crosses the cut: prefix head positions stay in [0, n-1]
            head, lo, hi = 0, 0, 0
            for m in word:
                if m == "L":
                    head += 1
                elif m == "R":
                    head -= 1
                lo, hi = min(lo, head), max(hi, head)
            if lo < 0 or hi > n - 1:
                ok = False
                log(f"FAIL part A (cut crossed): n={n} w={w}")
                break
            cnt += 1
        rows.append({"n": n, "lines": cnt, "N_rot": B - 1, "seconds": round(time.time() - t0, 1)})
        log(f"part A n={n}: {cnt} lines sorted, N_X = inv, N_rot = B_n - 1 = {B - 1}, {time.time() - t0:.1f} s")
    return ok, rows


def part_b(log):
    rows = []
    ok = True
    for n in range(4, 11, 2):
        path = os.path.join(ROOT, "data", "tables", f"dist_n{n}.bin")
        if not os.path.exists(path):
            continue
        circ = open(path, "rb").read()
        fact = factorials(n)
        m = n // 2
        tau = list(range(n))
        tau[m], tau[m + 1] = tau[m + 1], tau[m]
        tau = tuple(tau)
        d = circ[rank_perm(tau, fact)]
        # I(tau): minimum over all double cuts
        I = min(inversions(tuple((tau[(q + 1 + j) % n] - (q + 1 - c)) % n for j in range(n)))
                for q in range(n) for c in range(n))
        word = "L" * m + "X" + "R" * m
        sorts = apply_word(tau, word) == identity(n)
        good = (d == n + 1 and I == 1 and sorts and len(word) == n + 1)
        ok = ok and good
        rows.append({"n": n, "tau": list(tau), "d": d, "I": I, "word": word, "ok": good})
        log(f"part B n={n}: tau={tau} d={d} (n+1={n + 1}) I={I} word {word} sorts={sorts} -> {'ok' if good else 'FAIL'}")
    return ok, rows


def part_c(log):
    rows = []
    ok = True
    for n in range(4, 11):
        path = os.path.join(ROOT, "data", "runs", "line_profile", f"profile_n{n}.json")
        if not os.path.exists(path):
            continue
        r = json.load(open(path))
        maxI = max(row["I"] for row in r["by_I"])
        good = (r["max_d_minus_I"] == n * n // 4 and maxI == (n - 1) ** 2 // 4)
        ok = ok and good
        rows.append({"n": n, "max_d_minus_I": r["max_d_minus_I"], "floor_n2_4": n * n // 4,
                     "max_I": maxI, "floor_(n-1)2_4": (n - 1) ** 2 // 4, "max_d_minus_2I": r["max_d_minus_2I"]})
        log(f"part C n={n}: max(d - I) = {r['max_d_minus_I']} (floor(n^2/4) = {n * n // 4}), max I = {maxI} "
            f"(floor((n-1)^2/4) = {(n - 1) ** 2 // 4}), max(d - 2I) = {r['max_d_minus_2I']} -> {'ok' if good else 'FAIL'}")
    return ok, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amax", type=int, default=8)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    lines = []

    def log(msg):
        print(msg, flush=True)
        lines.append(msg)

    log(f"== {VERSION} {lc.CONSTRUCTION_VERSION} {CORE_VERSION} args={vars(args)}")
    oka, ra = part_a(args.amax, log)
    okb, rb = part_b(log)
    okc, rc = part_c(log)
    verdict = "PASS" if (oka and okb and okc) else "FAIL"
    log(f"verdict: {verdict}")
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump({"version": VERSION, "construction": lc.CONSTRUCTION_VERSION, "core": CORE_VERSION,
                   "args": vars(args), "part_a": ra, "part_b": rb, "part_c": rc, "verdict": verdict}, f, indent=1)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# check_C34 — линейная модель: универсальный маршрут, антиподальная транспозиция, профиль\n\n")
        f.write("Команда: `python3 checks/check_C34.py --amax %d`. Версии: %s, %s, %s.\n\n```text\n" %
                (args.amax, VERSION, lc.CONSTRUCTION_VERSION, CORE_VERSION))
        f.write("\n".join(lines) + "\n```\n")
    sys.exit(0 if verdict == "PASS" else 1)


if __name__ == "__main__":
    main()
