"""Checks for C28 (H10 on reflections, proved in docs/proofs/h10_reflections.md)
and finite checks of C27/H10 with a lean per-shift computation.

Core: oracle-1.0 (oracle/moves.py). Construction: strict_upper-1.0
(constructions/strict_upper.py); the carrier route is a variant of N1
(AGENTS.md rule 16), the lemmas of N1 are not transferred to it.

Everything is exact integer arithmetic.

Quantities per shift c (all of them are the ones of PROBLEM 5.2 / AGENTS rule 15):
    F_c, S_c, K_c (number of nontrivial cycles of f_c), carriers,
    H(c) = N1 external route, R_c = shortest walk 0 -> carriers -> c.
The value under test is `val(c) = 2F_c - S_c + R_c` (H10) and, where stated,
`valH(c) = 2F_c - S_c + H(c)` (the own route of N1, H9).

Modes:
    --reflections NMAX  closed forms of lemmas R1-R3 and the case analysis of the
                        theorem for every reflection pi(i) = h - i, 4 <= n <= NMAX,
                        every h, every c (exhaustive in h and c).
    --crosscheck NMAX   lean numbers == length of the actually built and executed
                        carrier-route word, all pi and all c, 4 <= n <= NMAX.
    --exhaustive N      all pi in S_N: max over pi of min_c val(c) - B_n
                        (checkpointed, resumable; writes data/runs/h10_reflections/).
    --census NMAX       tight inputs (min_c val(c) >= B_n) and the coverage of the
                        sufficient condition SC of the proof (Lemma G2), 4 <= n <= NMAX.

Usage:
    python3 checks/check_C28.py --reflections 40
    python3 checks/check_C28.py --crosscheck 7
    python3 checks/check_C28.py --census 8
    python3 checks/check_C28.py --exhaustive 10 [--chunk 200000]
"""

import argparse
import itertools
import json
import os
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
for sub in ("oracle", "constructions"):
    sys.path.insert(0, os.path.join(ROOT, sub))

from moves import delta, freely_reduce, CORE_VERSION  # noqa: E402
import strict_upper as su  # noqa: E402

CHECK_VERSION = "check_C28-1.0"
OUT_DIR = os.path.join(ROOT, "data", "runs", "h10_reflections")


def P(n):
    return (n * n) // 4


def B(n):
    return n * (n - 1) // 2


# ------------------------------------------------------------------ lean per-shift data

def shift_stats(pi, c):
    """(F_c, S_c, K_c, R_c, H(c), carriers) without building the word."""
    n = len(pi)
    f = [(pi[i] + c) % n for i in range(n)]
    F = S = K = 0
    carriers = []
    for cyc in su.cycles_of(f):
        cd = su.cycle_data(cyc, n)
        carriers.append(cyc[su.choose_carrier(cyc, cd, n)])
        F += cd["F"]
        S += cd["S"]
        K += 1
    R = len(su.carrier_route(carriers, c, n)[1])
    return F, S, K, R, su.H_of(c, n), carriers


def split_stats(pi, c):
    """(sum of M(C), E_c, N_X, N_rot) of the N1 word with the own route H(c);
    N_X and N_rot follow from the letter count of N1 2.4 (lemma 11)."""
    n = len(pi)
    f = [(pi[i] + c) % n for i in range(n)]
    sumM = E = F = K = 0
    for cyc in su.cycles_of(f):
        cd = su.cycle_data(cyc, n)
        sumM += cd["M"]
        E += cd["E"]
        F += cd["F"]
        K += 1
    return sumM, E, F - sumM - E + K, F - sumM + su.H_of(c, n)


def all_shifts(pi):
    return [shift_stats(pi, c) for c in range(len(pi))]


def min_excess(pi, route="carrier"):
    """min_c [2F_c - S_c + route_c] - B_n."""
    n = len(pi)
    best = None
    for (F, S, K, R, H, car) in all_shifts(pi):
        v = 2 * F - S + (R if route == "carrier" else H)
        best = v if best is None else min(best, v)
    return best - B(n)


# ------------------------------------------------------------------ mode: reflections

def reflection_forms(n, h, c):
    """Closed forms of lemmas R1-R3 for pi(i) = h - i: (F_c, S_c, K_c)."""
    if n % 2:
        K = (n - 1) // 2
        F = P(n)
    else:
        m = n // 2
        K = m - 1 if (h + c) % 2 == 0 else m
        if n % 4 == 0:
            F = P(n)
        else:
            F = P(n) - 1 if (c - h) % 2 == 0 else P(n) + 1
    return F, 3 * K, K


def theorem_shift(n, h):
    """The shift c chosen by the proof and the value of 2F_c - S_c + H(c) - B_n
    that the case analysis predicts."""
    if n % 2:
        return 1, 0
    if n % 4 == 0:
        return (1, -1) if h % 2 == 0 else (2, 0)
    return (1, 0) if h % 2 else (1, 1)


def reflection_min_H(n, h):
    """Closed form of min_c [2F_c - S_c + H(c)] - B_n on the reflection pi(i) = h - i
    (theorem 1 of docs/proofs/h10_reflections.md)."""
    if n % 2:
        return 0
    if n % 4 == 0:
        return -1 if h % 2 == 0 else 0
    return 1 if h % 2 == 0 else 0


def run_reflections(nmax):
    checked = 0
    for n in range(4, nmax + 1):
        for h in range(n):
            pi = tuple((h - i) % n for i in range(n))
            rows = all_shifts(pi)
            for c, (F, S, K, R, H, car) in enumerate(rows):
                pF, pS, pK = reflection_forms(n, h, c)
                assert (F, S, K) == (pF, pS, pK), ("R1-R3", n, h, c, (F, S, K), (pF, pS, pK))
                assert R <= H, ("R4", n, h, c, R, H)
                checked += 1
            cstar, pred = theorem_shift(n, h)
            F, S, K, R, H, car = rows[cstar]
            got = 2 * F - S + H - B(n)
            assert got == pred, ("case analysis", n, h, cstar, got, pred)
            assert got <= 1, ("theorem", n, h, got)
            # exact closed form of the minimum over c for the own route H(c)
            mnH = min(2 * F2 - S2 + H2 for (F2, S2, K2, R2, H2, _) in rows) - B(n)
            assert mnH == reflection_min_H(n, h), ("min_c H closed form", n, h, mnH)
            # the carrier route can only help
            assert min(2 * F2 - S2 + R2 for (F2, S2, K2, R2, H2, _) in rows) - B(n) <= 1, ("H10", n, h)
            # sharpness: n = 2 mod 4 and even h is the only case with +1
            mn = min(2 * F2 - S2 + R2 for (F2, S2, K2, R2, H2, _) in rows) - B(n)
            if mn == 1:
                assert n % 4 == 2 and h % 2 == 0, ("unexpected +1", n, h)
            if h == 1:
                # proposition 1 (C25 / H8): sigma_n, exact minimum B_n at the
                # stated shifts, with the geodesic split of C13
                mins = [c2 for c2, (F2, S2, K2, R2, H2, _) in enumerate(rows)
                        if 2 * F2 - S2 + H2 == B(n)]
                assert mnH == 0 and mins == sorted({0, 2, n - 2} if n % 4 == 0
                                                   else {1, n - 1}), ("C25 shifts", n, mins)
                for c2 in mins:
                    sumM, E, NX, NR = split_stats(pi, c2)
                    assert (NX, NR) == ((n - 1) ** 2 // 4, n * n // 4), ("C25 split", n, c2, NX, NR)
    print(f"reflections 4<=n<={nmax}: {checked} (h,c) pairs, lemmas R1-R4 and the "
          f"case analysis hold; min_c value <= B_n + 1 everywhere")
    return checked


# ------------------------------------------------------------------ mode: crosscheck

def run_crosscheck(nmax):
    cnt = 0
    for n in range(4, nmax + 1):
        for pi in itertools.permutations(range(n)):
            for c in range(n):
                F, S, K, R, H, car = shift_stats(pi, c)
                sw = su.shift_word(pi, c, "carrier")
                st = sw["stats"]
                assert (F, S, K) == (st["F_c"], st["S_c"], st["n_cycles"]), (pi, c)
                assert R == st["route_len"] and st["len"] == 2 * F - S + R, (pi, c)
                # lemma 11: letter count of the N1 word with its own route
                sumM, E, NX, NR = split_stats(pi, c)
                st1 = su.shift_word(pi, c, "n1")["stats"]
                assert (NX, NR) == (st1["N_X"], st1["N_rot"]), (pi, c, NX, NR)
                cnt += 1
    print(f"crosscheck 4<=n<={nmax}: {cnt} shifts, lean numbers == executed word length")
    return cnt


# ------------------------------------------------------------------ mode: census

def sufficient_condition(pi):
    """Lemma G2 (proof-usable bounds only): is there c with
    3K_c + 2(P - F_c) >= floor(n/2) - 1 + H(c)?"""
    n = len(pi)
    for (F, S, K, R, H, car) in all_shifts(pi):
        if 3 * K + 2 * (P(n) - F) >= n // 2 - 1 + H:
            return True
    return False


def is_reflection(pi):
    n = len(pi)
    hs = {(pi[i] + i) % n for i in range(n)}
    return hs.pop() if len(hs) == 1 else None


def run_census(nmax):
    out = {}
    for n in range(4, nmax + 1):
        hist = {}
        tight = []
        covered = 0
        total = 0
        for pi in itertools.permutations(range(n)):
            total += 1
            e = min_excess(pi)
            hist[e] = hist.get(e, 0) + 1
            if e >= 0:
                h = is_reflection(pi)
                tight.append({"pi": list(pi), "excess": e,
                              "reflection_h": h})
            if sufficient_condition(pi):
                covered += 1
        out[n] = {"total": total, "max_excess": max(hist), "hist": hist,
                  "tight": tight, "covered_by_SC": covered}
        print(f"n={n}: max min_c excess {max(hist)}, tight inputs {len(tight)} "
              f"({sum(1 for t in tight if t['reflection_h'] is not None)} reflections), "
              f"SC covers {covered}/{total} = {100.0*covered/total:.1f}%")
    return out


# ------------------------------------------------------------------ mode: exhaustive

def run_exhaustive(n, chunk, out_path):
    """All pi in S_n, checkpointed every `chunk` permutations. Resumable: reads
    out_path and continues from the stored index."""
    state = {"n": n, "done": 0, "hist": {}, "tight": [], "seconds": 0.0,
             "core": CORE_VERSION, "construction": su.CONSTRUCTION_VERSION,
             "check": CHECK_VERSION, "finished": False}
    if os.path.exists(out_path):
        with open(out_path) as fh:
            state = json.load(fh)
        state["hist"] = {int(k): v for k, v in state["hist"].items()}
        print(f"resume at {state['done']}")
    hist = {int(k): v for k, v in state["hist"].items()}
    it = itertools.permutations(range(n))
    for _ in range(state["done"]):
        next(it)
    t0 = time.time()
    done = state["done"]
    for pi in it:
        e = min_excess(pi)
        hist[e] = hist.get(e, 0) + 1
        if e >= 0:
            state["tight"].append({"pi": list(pi), "excess": e,
                                   "reflection_h": is_reflection(pi)})
        done += 1
        if done % chunk == 0:
            state.update({"done": done, "hist": hist,
                          "seconds": state["seconds"] + time.time() - t0})
            t0 = time.time()
            with open(out_path, "w") as fh:
                json.dump(state, fh)
            print(f"  {done} done, max excess {max(hist)}, {state['seconds']:.0f} s", flush=True)
    state.update({"done": done, "hist": hist, "finished": True,
                  "seconds": state["seconds"] + time.time() - t0,
                  "max_excess": max(hist)})
    with open(out_path, "w") as fh:
        json.dump(state, fh)
    print(f"n={n}: {done} permutations, max min_c excess {max(hist)}, "
          f"tight inputs {len(state['tight'])}, {state['seconds']:.0f} s")
    return state


# ------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reflections", type=int, default=0)
    ap.add_argument("--crosscheck", type=int, default=0)
    ap.add_argument("--census", type=int, default=0)
    ap.add_argument("--exhaustive", type=int, default=0)
    ap.add_argument("--chunk", type=int, default=200000)
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)
    print(CHECK_VERSION, CORE_VERSION, su.CONSTRUCTION_VERSION)
    if args.reflections:
        run_reflections(args.reflections)
    if args.crosscheck:
        run_crosscheck(args.crosscheck)
    if args.census:
        res = run_census(args.census)
        with open(os.path.join(OUT_DIR, "census.json"), "w") as fh:
            json.dump(res, fh, indent=1)
    if args.exhaustive:
        run_exhaustive(args.exhaustive, args.chunk,
                       os.path.join(OUT_DIR, f"exhaustive_n{args.exhaustive}.json"))


if __name__ == "__main__":
    main()
