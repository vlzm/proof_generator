"""Checker for C37 (session 9): exact one-step recursion for value-rotation inversions.

Claim (proved in docs/notes/h13_line_model.md Sec. 7, Lemma D): for any sequence
v = (v_0, ..., v_{n-1}) that is a permutation of {0, ..., n-1}, and any b in
Z_n, define u^(b)_j = (v_j - b) mod n and inv(b) = #{(j,k): j<k, u^(b)_j >
u^(b)_k}.  Then

    inv(b+1 mod n) = inv(b) + (n - 1 - 2 * pos(v, b))

where pos(v, b) is the index j with v_j = b (i.e. pos(v, b) = v^{-1}(b)).

Consequence used in the note: writing this as a prefix-sum identity over
P = v^{-1} recovers inv(v) = inv(v^{-1}) (standard fact, re-derived here as a
corollary, not assumed) and reformulates min_b inv(b) as a prefix-sum
extremal problem on P.  The note also records (Sec. 7.2) that restricting to
this single value-rotation (i.e. fixing the position order, freedom in b
only) is NOT enough to prove H13-I: it fails exactly at the two inputs
already on record as counterexamples to the weaker "locally optimal cut"
condition (Lemma C in the note): n = 7, v = (0,5,4,3,2,1,6) and n = 8,
v = (0,6,5,3,4,2,1,7), where min_b inv(b) exceeds floor((n-1)^2/4) by 1.  This
checker re-verifies both facts: the recursion (exhaustively for small n,
randomly for larger n) and the two recorded failures of the value-rotation-only
bound, plus inv(v) = inv(v^{-1}) exhaustively.

H13-I itself (existence of a *joint* position+value rotation with
inv <= floor((n-1)^2/4) for every pi) is NOT decided by this checker: it
remains CONJECTURED / VERIFIED only for 4 <= n <= 10 (C33).

Usage: python3 checks/check_C37.py [--emax 8] [--rand_nmax 200] [--rand_trials 2000] [--seed 0]
Output: data/runs/check_C37/report.json, report.md.
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
from moves import CORE_VERSION  # noqa: E402  (only for the version stamp; no moves used here)

VERSION = "check_C37-1.0"
OUT = os.path.join(ROOT, "data", "runs", "check_C37")


def inv_count(u):
    """O(n log n) merge-sort inversion count (independent of the recursion under test)."""
    if len(u) <= 1:
        return 0
    a = list(u)

    def sort_count(lo, hi):
        if hi - lo <= 1:
            return 0
        mid = (lo + hi) // 2
        c = sort_count(lo, mid) + sort_count(mid, hi)
        merged = []
        i, j = lo, mid
        while i < mid and j < hi:
            if a[i] <= a[j]:
                merged.append(a[i])
                i += 1
            else:
                merged.append(a[j])
                j += 1
                c += mid - i
        merged.extend(a[i:mid])
        merged.extend(a[j:hi])
        a[lo:hi] = merged
        return c

    return sort_count(0, len(a))


def inv_all_b(v):
    """Direct computation of inv(b) for all b, for cross-checking the recursion."""
    n = len(v)
    return [inv_count([(x - b) % n for x in v]) for b in range(n)]


def recursion_check(v, log, ctx):
    """Verify inv(b+1) = inv(b) + (n-1-2*pos(v,b)) for all b, against direct computation."""
    n = len(v)
    direct = inv_all_b(v)
    pos = [0] * n
    for j, val in enumerate(v):
        pos[val] = j
    ok = True
    for b in range(n):
        predicted = direct[b] + (n - 1 - 2 * pos[b])
        actual = direct[(b + 1) % n]
        if predicted != actual:
            ok = False
            log(f"FAIL recursion: {ctx} v={v} b={b} pos(v,b)={pos[b]} "
                f"predicted inv(b+1)={predicted} actual={actual}")
    return ok


def part_recursion_exhaustive(emax, log):
    rows = []
    ok = True
    for n in range(2, emax + 1):
        t0 = time.time()
        cnt = 0
        for v in itertools.permutations(range(n)):
            if not recursion_check(list(v), log, f"exhaustive n={n}"):
                ok = False
            cnt += 1
        rows.append({"n": n, "perms": cnt, "seconds": round(time.time() - t0, 2)})
        log(f"recursion exhaustive n={n}: {cnt} permutations, {time.time() - t0:.2f} s")
    return ok, rows


def part_recursion_random(rand_nmax, trials, seed, log):
    rng = random.Random(seed)
    rows = []
    ok = True
    ns = sorted(set([9, 10, 20, 50] + ([rand_nmax] if rand_nmax not in (9, 10, 20, 50) else [])))
    ns = [n for n in ns if n <= rand_nmax]
    for n in ns:
        t0 = time.time()
        checked = 0
        for _ in range(trials):
            v = list(range(n))
            rng.shuffle(v)
            if not recursion_check(v, log, f"random n={n}"):
                ok = False
            checked += 1
        rows.append({"n": n, "trials": checked, "seconds": round(time.time() - t0, 2)})
        log(f"recursion random n={n}: {checked} trials, {time.time() - t0:.2f} s")
    return ok, rows


def part_inv_inverse_exhaustive(emax, log):
    """inv(v) == inv(v^{-1}) for every permutation, 2 <= n <= emax."""
    ok = True
    rows = []
    for n in range(2, emax + 1):
        cnt = 0
        for v in itertools.permutations(range(n)):
            vinv = [0] * n
            for j, val in enumerate(v):
                vinv[val] = j
            if inv_count(v) != inv_count(vinv):
                ok = False
                log(f"FAIL inv(v)=inv(v^-1): n={n} v={v}")
            cnt += 1
        rows.append({"n": n, "perms": cnt})
        log(f"inv(v)=inv(v^-1) n={n}: {cnt} permutations checked")
    return ok, rows


KNOWN_FAILURES = {
    7: (0, 5, 4, 3, 2, 1, 6),
    8: (0, 6, 5, 3, 4, 2, 1, 7),
}


def part_value_rotation_only_insufficient(log):
    """Re-verify the two recorded counterexamples to the value-rotation-only bound
    min_b inv(b) <= floor((n-1)^2/4) (position order fixed, only b free): both
    exceed the target by exactly 1, matching docs/notes/h13_line_model.md Sec. 4
    (Lemma C) and the session-9 addendum (Sec. 7.2)."""
    ok = True
    rows = []
    for n, v in KNOWN_FAILURES.items():
        target = (n - 1) ** 2 // 4
        m = min(inv_all_b(list(v)))
        good = (m == target + 1)
        ok = ok and good
        rows.append({"n": n, "v": list(v), "min_b_inv": m, "target": target, "ok": good})
        log(f"value-rotation-only n={n} v={v}: min_b inv = {m}, floor((n-1)^2/4) = {target}, "
            f"exceeds by {m - target} -> {'ok (matches recorded failure)' if good else 'FAIL (does not match record)'}")
    return ok, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--emax", type=int, default=8, help="exhaustive upper bound on n")
    ap.add_argument("--rand_nmax", type=int, default=200)
    ap.add_argument("--rand_trials", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    lines = []

    def log(msg):
        print(msg, flush=True)
        lines.append(msg)

    log(f"== {VERSION} {CORE_VERSION} args={vars(args)}")
    ok1, r1 = part_recursion_exhaustive(args.emax, log)
    ok2, r2 = part_recursion_random(args.rand_nmax, args.rand_trials, args.seed, log)
    ok3, r3 = part_inv_inverse_exhaustive(args.emax, log)
    ok4, r4 = part_value_rotation_only_insufficient(log)
    verdict = "PASS" if (ok1 and ok2 and ok3 and ok4) else "FAIL"
    log(f"verdict: {verdict}")

    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump({"version": VERSION, "core": CORE_VERSION, "args": vars(args),
                   "recursion_exhaustive": r1, "recursion_random": r2,
                   "inv_inverse_exhaustive": r3, "value_rotation_only": r4,
                   "verdict": verdict}, f, indent=1)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# check_C37 — рекуррента для инверсий одного вращения значений (H13-I, сессия 9)\n\n")
        f.write("Команда: `python3 checks/check_C37.py --emax %d --rand_nmax %d --rand_trials %d --seed %d`. "
                "Версии: %s, %s.\n\n```text\n" %
                (args.emax, args.rand_nmax, args.rand_trials, args.seed, VERSION, CORE_VERSION))
        f.write("\n".join(lines) + "\n```\n")
    sys.exit(0 if verdict == "PASS" else 1)


if __name__ == "__main__":
    main()
