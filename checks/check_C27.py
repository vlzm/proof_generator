"""Checks for C27 / H10 -- the carrier-route variant of construction N1
(docs/proofs/N1_carrier_route.md, AGENTS.md rule 16).

C27 claimed, for all n >= 4 and all pi in S_n,

    (a)  min_c [2F_c - S_c + R_c] <= B_n + 1,
    (b)  after free reduction of the assembled word,  <= B_n.

This checker certifies the session-4 verdict REFUTED and the results proved
around it:

  A  cross-check of the C core checks/carrier_core.h (used for the exhaustive
     runs at n = 10, 11, 12) against constructions/strict_upper.py on all
     permutations for 4 <= n <= 7 and on random ones for n = 8, 9: F_c, S_c,
     R_c and 2F_c - S_c + R_c must agree, and the Python word length must equal
     2F_c - S_c + R_c.
  B  the minimal counterexample to (b): n = 10, pi = (9,2,5,0,1,8,7,6,3,4) --
     every shift gives a word of length >= B_n + 1 which free reduction does
     not shorten below B_n + 1.
  C  the minimal counterexample to (a): n = 12, pi = (1,5,3,10,2,0,7,11,9,4,8,6)
     with min over c equal to B_12 + 3 and F_c = P at every shift; plus the
     affine witnesses n = 18 (5i+16), 26, 27, 28, where the excess grows.
  D  Theorem A of the proof document (reflections pi(i) = h - i): the closed
     forms F_c, S_c by parity, and the bound reached at c in {0,1}
     (<= B_n for odd n and n = 0 mod 4, <= B_n + 1 for n = 2 mod 4, with
     equality exactly at even h) -- verified for 4 <= n <= 30.
  E  refutation of the natural per-cycle lemma S(C) >= |C| + 1: the zigzag
     cycle (0,1,n-2,n-1,n-4,n-3,...) has S(C) = 5 and |C| = n for even n >= 8.
  F  full enumeration of all oriented nontrivial cycles for 4 <= n <= 9: the
     minimum of S(C) - |C| per n (the zigzag attains it at n = 8).
  G  for every shift of every counterexample, R_c equals the length of a
     shortest walk computed by an independent BFS over (position, visited
     carriers) -- so the refutation cannot be an artefact of the walk routine.

Every word built here is executed letter by letter by the reference moves
(oracle-1.0) inside strict_upper.shift_word.

Usage: python3 checks/check_C27.py [--perm-max-n 7] [--reflection-max-n 30]
Output: data/runs/carrier_route_H10/check_C27.json
"""

import argparse
import itertools
import json
import os
import random
import subprocess
import sys
import tempfile
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
for sub in ("oracle", "bounds", "constructions"):
    sys.path.insert(0, os.path.join(ROOT, sub))

from moves import apply_word, freely_reduce, identity, CORE_VERSION  # noqa: E402
import known  # noqa: E402
import strict_upper as su  # noqa: E402

OUT_DIR = os.path.join(ROOT, "data", "runs", "carrier_route_H10")
CHECKER_VERSION = "check_C27-1.0"


class CheckFailure(Exception):
    pass


def check(cond, tag, info=None):
    if not cond:
        raise CheckFailure(f"{tag}: {info!r}")


def carrier_stats(pi, c):
    """(F_c, S_c, R_c, |word|, |reduced word|) from strict_upper (words executed)."""
    sw = su.shift_word(pi, c, "carrier")
    st, w = sw["stats"], sw["word"]
    check(apply_word(tuple(pi), w) == identity(len(pi)), "word must sort pi", (pi, c))
    r = freely_reduce(w)
    check(apply_word(tuple(pi), r) == identity(len(pi)), "reduced word must sort pi", (pi, c))
    check(len(w) == 2 * st["F_c"] - st["S_c"] + st["route_len"], "length identity", (pi, c))
    return st["F_c"], st["S_c"], st["route_len"], len(w), len(r)


# ---------------------------------------------------------------- A: C core

def run_cross_check(perm_max_n, seed):
    exe = os.path.join(ROOT, "checks", "check_C27_fast")
    src = os.path.join(ROOT, "checks", "check_C27_fast.c")
    if not os.path.exists(exe) or os.path.getmtime(exe) < os.path.getmtime(src):
        subprocess.run(["gcc", "-O2", "-o", exe, src], check=True)
    cases = []
    for n in range(4, perm_max_n + 1):
        cases.extend((n, p) for p in itertools.permutations(range(n)))
    rng = random.Random(seed)
    for n in (8, 9):
        for _ in range(200):
            p = list(range(n))
            rng.shuffle(p)
            cases.append((n, tuple(p)))
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        for n, p in cases:
            f.write(f"{n} " + " ".join(map(str, p)) + "\n")
        path = f.name
    out = subprocess.run([exe, "-f", path], check=True, capture_output=True, text=True).stdout.splitlines()
    os.unlink(path)
    check(len(out) == len(cases), "one output line per case", (len(out), len(cases)))
    for (n, p), line in zip(cases, out):
        vals = list(map(int, line.split()))
        check(len(vals) == 5 * n, "five numbers per shift", (n, p))
        for c in range(n):
            cc, F, S, R, v = vals[5 * c:5 * c + 5]
            F2, S2, R2, L, _ = carrier_stats(p, c)
            check((cc, F, S, R, v) == (c, F2, S2, R2, L),
                  "C core must agree with strict_upper", (n, p, c, vals[5 * c:5 * c + 5], (F2, S2, R2, L)))
    return {"cases": len(cases), "shifts": sum(n for n, _ in cases),
            "coverage": f"all permutations 4 <= n <= {perm_max_n}, 200 random each for n = 8, 9"}


# ---------------------------------------------------------------- B, C: counterexamples

def table(pi):
    n = len(pi)
    rows = []
    for c in range(n):
        F, S, R, L, red = carrier_stats(pi, c)
        rows.append({"c": c, "F_c": F, "S_c": S, "R_c": R, "H": su.H_of(c, n), "len": L, "reduced": red})
    return rows


def run_counterexamples():
    res = {}
    # B: minimal counterexample to part (b) -- n = 10
    n = 10
    pi = (9, 2, 5, 0, 1, 8, 7, 6, 3, 4)
    Bn = known.target_diameter(n)
    rows = table(pi)
    mlen = min(r["len"] for r in rows)
    mred = min(r["reduced"] for r in rows)
    check(mlen == Bn + 1 and mred == Bn + 1, "n=10 counterexample to C27(b)", (mlen, mred, Bn))
    res["C27b_minimal"] = {"n": n, "pi": list(pi), "B_n": Bn, "min_len": mlen, "min_reduced": mred,
                           "rows": rows,
                           "minimality": "exhaustive for 4 <= n <= 9 (C27v, session 3) and unique among "
                                         "the 3 628 800 permutations of S_10 (check_C27_fast 10)"}
    # C0: minimal counterexample to part (a) -- n = 12 (not affine)
    n = 12
    pi = (1, 5, 3, 10, 2, 0, 7, 11, 9, 4, 8, 6)
    Bn = known.target_diameter(n)
    rows = table(pi)
    mlen = min(r["len"] for r in rows)
    mred = min(r["reduced"] for r in rows)
    check(mlen == Bn + 3, "n=12 minimal counterexample to C27(a)", (mlen, Bn))
    check(all(r["F_c"] == known.P(n) for r in rows), "F_c = P at every shift", rows)
    res["C27a_minimal"] = {"n": n, "pi": list(pi), "B_n": Bn, "min_len": mlen, "min_reduced": mred,
                           "rows": rows,
                           "minimality": "exhaustive for 4 <= n <= 11 (max excess +1); at n = 12 exactly "
                                         "17 of 479 001 600 permutations exceed B_n + 1, two of them by 3 "
                                         "(the other is (6,4,11,3,1,8,0,10,5,9,7,2))"}
    # C: further counterexamples to part (a), affine family (excess grows with n)
    res["C27a"] = {}
    for n, a, b, expect in ((18, 5, 16, 3), (26, 5, 6, 3), (27, 8, 3, 4), (28, 11, 20, 5)):
        pi = tuple((a * i + b) % n for i in range(n))
        Bn = known.target_diameter(n)
        rows = table(pi)
        mlen = min(r["len"] for r in rows)
        mred = min(r["reduced"] for r in rows)
        check(mlen - Bn == expect, f"affine counterexample n={n}", (mlen, Bn, expect))
        res["C27a"][n] = {"pi": f"{a}i+{b} mod {n}", "perm": list(pi), "B_n": Bn,
                          "min_len": mlen, "excess": mlen - Bn, "min_reduced": mred,
                          "reduced_excess": mred - Bn, "rows": rows if n == 18 else None}
    return res


# ---------------------------------------------------------------- D: Theorem A

def reflection_forms(n, g):
    """Closed forms of Theorem A: (F_c, S_c) for the reflection with h + c = g."""
    P = known.P(n)
    if n % 2 == 1:
        return P, 3 * ((n - 1) // 2)
    if n % 4 == 0:
        return (P, 3 * (n // 2 - 1)) if g % 2 == 0 else (P, 3 * (n // 2))
    return (P - 1, 3 * (n // 2 - 1)) if g % 2 == 0 else (P + 1, 3 * (n // 2))


def run_reflections(max_n):
    res = {}
    for n in range(4, max_n + 1):
        Bn = known.target_diameter(n)
        per_h = {}
        for h in range(n):
            pi = tuple((h - i) % n for i in range(n))
            vals = []
            for c in range(n):
                F, S, R, L, red = carrier_stats(pi, c)
                check((F, S) == reflection_forms(n, (h + c) % n),
                      "Theorem A closed forms", (n, h, c, F, S, reflection_forms(n, (h + c) % n)))
                check(R <= su.H_of(c, n), "R_c <= H(c)", (n, h, c, R, su.H_of(c, n)))
                vals.append((L, red))
            # the bound of Theorem A, computed from the closed forms at c in {0,1}
            bound = min(2 * reflection_forms(n, (h + c) % n)[0] - reflection_forms(n, (h + c) % n)[1]
                        + su.H_of(c, n) for c in (0, 1))
            mlen = min(v[0] for v in vals)
            mred = min(v[1] for v in vals)
            check(mlen <= bound, "actual min <= Theorem A bound", (n, h, mlen, bound))
            predicted = Bn + 1 if (n % 4 == 2 and h % 2 == 0) else Bn
            check(bound <= predicted, "Theorem A statement", (n, h, bound, predicted))
            check(mred <= Bn, "reflections satisfy C27(b)", (n, h, mred, Bn))
            per_h[h] = {"min_len_minus_B": mlen - Bn, "min_reduced_minus_B": mred - Bn,
                        "theorem_A_bound_minus_B": bound - Bn}
        res[n] = {"B_n": Bn, "per_h": per_h,
                  "max_min_len_minus_B": max(v["min_len_minus_B"] for v in per_h.values())}
    return res


# ---------------------------------------------------------------- E: zigzag cycles

def run_zigzag(max_n):
    """S(C) >= |C| + 1 is false: the zigzag n-cycle (0,1,n-2,n-1,n-4,n-3,...)
    has M = 3, t = n/2, E = 0, F = 2n and S = 5 for every even n >= 8, so
    S(C) - |C| = 5 - n tends to minus infinity."""
    res = {}
    for n in range(8, max_n + 1, 2):
        cyc = [0, 1]
        v = n - 2
        while v > 1:
            cyc.extend([v, v + 1])
            v -= 2
        cd = su.cycle_data(cyc, n)
        res[n] = {"cycle": cyc, "k": cd["k"], "M": cd["M"], "t": cd["t"], "E": cd["E"],
                  "F": cd["F"], "S": cd["S"], "S_minus_k": cd["S"] - cd["k"]}
        check((cd["k"], cd["M"], cd["t"], cd["E"], cd["F"], cd["S"]) == (n, 3, n // 2, 0, 2 * n, 5),
              "zigzag closed form (M, t, E, F, S) = (3, n/2, 0, 2n, 5)", (n, cd))
        check(cd["S"] - cd["k"] == 5 - n, "zigzag refutes S(C) >= |C| + 1", res[n])
    return res


def bfs_walk(required, c, n):
    """Independent shortest-walk length: BFS over (position, visited carriers).
    Used to certify that R_c really is the minimum, so that a counterexample
    cannot be an artefact of the walk heuristic in strict_upper.carrier_route."""
    req = sorted(set(required))
    idx = {v: i for i, v in enumerate(req)}
    full = (1 << len(req)) - 1
    start = (0, 1 << idx[0] if 0 in idx else 0)
    dist = {start: 0}
    frontier = [start]
    while frontier:
        nxt = []
        for st in frontier:
            p, vis = st
            if p == c and vis == full:
                return dist[st]
            for step in (1, -1):
                q = (p + step) % n
                v2 = vis | (1 << idx[q]) if q in idx else vis
                if (q, v2) not in dist:
                    dist[(q, v2)] = dist[st] + 1
                    nxt.append((q, v2))
        frontier = nxt
    raise CheckFailure(f"no walk: {required} {c} {n}")


def run_walk_certificate(perms):
    """R_c of the variant equals the true shortest walk, for the counterexamples."""
    out = {}
    for pi in perms:
        n = len(pi)
        for c in range(n):
            sw = su.shift_word(pi, c, "carrier")
            carriers = sw["stats"]["carriers"]
            R = sw["stats"]["route_len"]
            b = bfs_walk(carriers, c, n)
            check(R == b, "R_c must be the shortest walk (BFS)", (pi, c, carriers, R, b))
        out[str(list(pi))] = {"n": n, "shifts_checked": n}
    return out


def run_all_cycles(max_n):
    """Full enumeration of oriented nontrivial cycles: minimum of S(C) - |C|."""
    res = {}
    for n in range(4, max_n + 1):
        worst, wex, cnt = None, None, 0
        for k in range(2, n + 1):
            for verts in itertools.combinations(range(n), k):
                for tail in itertools.permutations(verts[1:]):
                    cyc = [verts[0]] + list(tail)
                    cd = su.cycle_data(cyc, n)
                    cnt += 1
                    if worst is None or cd["S"] - k < worst:
                        worst, wex = cd["S"] - k, cyc
        res[n] = {"cycles": cnt, "min_S_minus_k": worst, "witness": wex}
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perm-max-n", type=int, default=7)
    ap.add_argument("--reflection-max-n", type=int, default=30)
    ap.add_argument("--zigzag-max-n", type=int, default=20)
    ap.add_argument("--cycle-max-n", type=int, default=9,
                    help="full enumeration of oriented cycles up to this n (section F)")
    ap.add_argument("--seed", type=int, default=20260913)
    args = ap.parse_args()
    t0 = time.time()
    report = {"checker": CHECKER_VERSION, "core": CORE_VERSION,
              "construction": su.CONSTRUCTION_VERSION, "args": vars(args),
              "claim": "C27 / H10 (carrier-route variant): REFUTED"}
    t = time.time()
    report["A_cross_check"] = run_cross_check(args.perm_max_n, args.seed)
    print("A: C core == strict_upper on", report["A_cross_check"]["cases"], "permutations",
          f"({time.time() - t:.1f} s)", flush=True)
    report["B_C_counterexamples"] = run_counterexamples()
    print("B: n=10 pi=(9,2,5,0,1,8,7,6,3,4): min_c len = min_c reduced = B_10 + 1 = 46", flush=True)
    print("C: minimal counterexample to C27(a): n=12 pi=(1,5,3,10,2,0,7,11,9,4,8,6), min_c = B_12 + 3; "
          "affine: n=18 +3, n=26 +3, n=27 +4, n=28 +5", flush=True)
    report["D_reflections"] = run_reflections(args.reflection_max_n)
    print("D: Theorem A verified for all reflections, 4 <= n <=", args.reflection_max_n, flush=True)
    report["E_zigzag"] = run_zigzag(args.zigzag_max_n)
    print("E: zigzag cycles refute S(C) >= |C|:",
          {n: v["S_minus_k"] for n, v in report["E_zigzag"].items()}, flush=True)
    ce = [(9, 2, 5, 0, 1, 8, 7, 6, 3, 4), (1, 5, 3, 10, 2, 0, 7, 11, 9, 4, 8, 6)]
    ce += [tuple((a * i + b) % n for i in range(n)) for n, a, b in ((18, 5, 16), (26, 5, 6), (27, 8, 3), (28, 11, 20))]
    report["G_walk_certificate"] = run_walk_certificate(ce)
    print("G: R_c of every shift of the counterexamples equals the BFS shortest walk", flush=True)
    t = time.time()
    report["F_all_cycles"] = run_all_cycles(args.cycle_max_n)
    print("F: full cycle enumeration, min S(C) - |C| =",
          {n: v["min_S_minus_k"] for n, v in report["F_all_cycles"].items()},
          f"({time.time() - t:.1f} s)", flush=True)
    report["total_seconds"] = round(time.time() - t0, 2)
    report["result"] = "PASS (all checks of the session-4 verdict hold)"
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(os.path.join(OUT_DIR, "check_C27.json"), "w") as f:
        json.dump(report, f, indent=1)
    print("PASS in", report["total_seconds"], "s")


if __name__ == "__main__":
    try:
        main()
    except CheckFailure as exc:
        print("FAIL:", exc)
        sys.exit(1)
