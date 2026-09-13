"""Finite checks for C32 (docs/notes/h12_verdict.md): the "horizon" head policy
(constructions/candidate_sweep.py, sweep-1.3, mode "horizon", horizon=3,
weight=0.5, thr in {1, 2}) on reflections pi_h(i) = h - i mod n, restricted to
shift c in {0, 1}.

This is NOT a proof: it is the numeric evidence behind the C32 row and behind
the h12_verdict.md decision not to attempt a formal proof (AGENTS.md rule 5 -
finite search does not certify the general claim). Checked facts, for every
h in Z_n and 4 <= n <= --nmax:

  1. Exactly one of c in {0, 1} avoids deadlock at n even (both do at n odd),
     and the working c always matches the C28 lemma-4 parity rule: h' = h + c
     is odd when n = 0 (mod 4), even when n = 2 (mod 4).
  2. Every X in the resulting word is either a "double" swap (both elements
     move closer to target, Phi drops by 2) or a "single"/bystander swap
     (Phi unchanged: one element is carried past a position currently at its
     own target, rule 15 "single"); no swap changes Phi by +-1 or +2 (the
     coded gain heuristic never fires a net-harmful swap here).
  3. double == F_c / 2 with F_c from C28 lemma 2 at the working shift's h'
     parity; single == floor((n - 2)^2 / 8), the SAME value for every h at
     fixed n (an invariance under head-start rotation that is not proved -
     see h12_verdict.md); N_X == double + single == floor((n - 1)^2 / 4).

Word replay uses only the (word, c) pair returned by candidate_sweep and a
direct recomputation of rem/Phi (oracle definitions, PROBLEM section 2 /
rule 16) - it does not read candidate_sweep internals.

Usage: python3 checks/check_C32.py [--nmax 30] [--large 50,60,80] [--large-samples 8]
Runtime grows roughly as O(n^4) per h (bounded-lookahead search rescanned
per decision point), so nmax and --large are kept modest by default -
nmax = 30 exhaustive is already ~140 s, n = 80 sampled ~200-300 s; raise
only with a checkpoint (AGENTS.md rule 10).
Output: data/runs/check_C32/report.json
"""

import argparse
import json
import os
import random
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
for sub in ("oracle", "constructions"):
    sys.path.insert(0, os.path.join(ROOT, sub))

from moves import apply_word, identity, CORE_VERSION  # noqa: E402
import candidate_sweep as cs  # noqa: E402

CFG = dict(horizon=3, weight=0.5, thrs=(1, 2))


def refl(n, h):
    return tuple((h - i) % n for i in range(n))


def target_NX(n):
    return ((n - 1) * (n - 1)) // 4


def single_formula(n):
    return ((n - 2) * (n - 2)) // 8


def Fc_formula(n, hp):
    """C28 lemma 2, evaluated at h' = hp."""
    if n % 2 == 1:
        return (n * n - 1) // 4
    if hp % 2 == 0:
        return (n * n) // 4 if n % 4 == 0 else (n * n - 4) // 4
    return (n * n) // 4 if n % 4 == 0 else (n * n + 4) // 4


def good_parity(n, hp):
    """C28 lemma 4 deadlock-avoiding parity of h' = h + c."""
    if n % 2 == 1:
        return True
    if n % 4 == 0:
        return hp % 2 == 1
    return hp % 2 == 0


def replay_breakdown(pi, c, word):
    """Replay word with an independent rem/Phi recomputation (not candidate_sweep
    internals). Returns (double, single, other, final_circle)."""
    n = len(pi)
    target = [(e + c) % n for e in range(n)]
    circle = list(pi)

    def rem(pos):
        d = (target[circle[pos]] - pos) % n
        if 2 * d > n:
            d -= n
        return d

    head = 0
    double = single = other = 0
    for ch in word:
        if ch == "X":
            a, b = head, (head + 1) % n
            before = abs(rem(a)) + abs(rem(b))
            circle[a], circle[b] = circle[b], circle[a]
            after = abs(rem(a)) + abs(rem(b))
            d = after - before
            if d == -2:
                double += 1
            elif d == 0:
                single += 1
            else:
                other += 1
        elif ch == "L":
            head = (head + 1) % n
        elif ch == "R":
            head = (head - 1) % n
        else:
            raise ValueError(ch)
    return double, single, other, circle


def check_one(n, h):
    pi = refl(n, h)
    # which of c in {0, 1} avoid deadlock, and does the working one match lemma 4?
    works = {}
    for c in (0, 1):
        ok = False
        for dirn in (1, -1):
            for thr in CFG["thrs"]:
                try:
                    cs.sweep_word(pi, c, dirn, thr, mode="horizon", horizon=CFG["horizon"],
                                   weight=CFG["weight"], fallback=False)
                    ok = True
                except cs.ConstructionError:
                    continue
        works[c] = ok
    if n % 2 == 1:
        assert works[0] and works[1], ("odd n: both c should work", n, h)
    else:
        assert works[0] != works[1], ("even n: exactly one of c in {0,1} should work", n, h)
        good_c = 0 if works[0] else 1
        hp = (h + good_c) % n
        assert good_parity(n, hp), ("working c parity mismatches C28 lemma 4", n, h, good_c, hp)

    r = cs.best_sweep(pi, thrs=CFG["thrs"], modes=("horizon",), horizon=CFG["horizon"],
                       weight=CFG["weight"], cs=[0, 1], fallback=False)
    c = r["c"]
    hp = (h + c) % n
    double, single, other, final_circle = replay_breakdown(pi, c, r["word"])
    expected_final = [(p - c) % n for p in range(n)]
    assert final_circle == expected_final, ("replay does not sort", n, h)
    assert apply_word(tuple(pi), r["word"]) == identity(n), ("oracle: word does not sort", n, h)
    assert other == 0, ("swap with |dPhi| != 0, 2", n, h, other)
    assert double + single == r["N_X"], ("replay N_X mismatch", n, h, double, single, r["N_X"])
    exp_double = Fc_formula(n, hp) // 2
    exp_single = single_formula(n)
    assert double == exp_double, ("double != F_c/2", n, h, double, exp_double)
    assert single == exp_single, ("single != floor((n-2)^2/8)", n, h, single, exp_single)
    assert r["N_X"] == target_NX(n), ("N_X != floor((n-1)^2/4)", n, h, r["N_X"], target_NX(n))
    return double, single


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=30)
    ap.add_argument("--large", type=str, default="50,60")
    ap.add_argument("--large-samples", type=int, default=8)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    t0 = time.time()
    checked = 0
    for n in range(4, a.nmax + 1):
        for h in range(n):
            check_one(n, h)
            checked += 1
    exhaustive_s = time.time() - t0
    print(f"exhaustive PASS: 4<=n<={a.nmax}, all h, {checked} instances, {exhaustive_s:.1f} s")

    rng = random.Random(a.seed)
    large_rows = []
    for n in [int(x) for x in a.large.split(",") if x]:
        t1 = time.time()
        hs = sorted(set([0, 1, 2, 3, n // 4, n // 4 + 1, n // 2, n // 2 + 1, n - 2, n - 1] +
                         rng.sample(range(n), min(a.large_samples, n))))
        for h in hs:
            check_one(n, h)
        row = {"n": n, "h_sampled": hs, "seconds": round(time.time() - t1, 1)}
        large_rows.append(row)
        print(f"n={n}: SAMPLED {len(hs)} h, PASS, {row['seconds']} s")

    report = {
        "claim": "C32", "core": CORE_VERSION, "candidate": cs.CANDIDATE_VERSION,
        "cfg": CFG,
        "exhaustive_range": f"4<=n<={a.nmax}, all h in Z_n",
        "exhaustive_instances": checked, "exhaustive_seconds": round(exhaustive_s, 1),
        "large_sampled": large_rows,
        "result": "PASS (numeric evidence only, not a proof - see docs/notes/h12_verdict.md)",
    }
    out = os.path.join(ROOT, "data", "runs", "check_C32")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "report.json"), "w") as f:
        json.dump(report, f, indent=1)
    print(f"written {os.path.join(out, 'report.json')}")


if __name__ == "__main__":
    main()
