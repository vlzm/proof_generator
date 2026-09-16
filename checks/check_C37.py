"""Checker for C37 (double-cut identities), C38 (row-sum strengthening of H13-I)
and C39 (refuted proof strategies for H13-I), session 9.

Model.  For a permutation pi of Z_n and a double cut (a, b) in Z_n x Z_n let
    w_j(a, b) = (pi(a + j) - b) mod n,  j = 0..n-1,
    inv(a, b) = inversions of w(a, b),  I(pi) = min over the n^2 cuts,
    B_n = n(n-1)/2 = C(n,2),  S(a, b) = B_n - 2 inv(a, b).
H13-I is I(pi) <= floor((n-1)^2/4), equivalently max_{a,b} S(a,b) >= floor(n/2).

Part 1 (C37, docs/proofs/C37_torus_cut_identities.md): all lemmas, exhaustively.
  L1  S = concordant - discordant, and the threshold floor(n/2);
  L2  concordant pairs = agreements of the two linearizations (cyclic orders);
      reflections give exactly floor((n-1)^2/4);
  L3  inv(a+1,b) - inv(a,b) = (n-1) - 2 w_0,  inv(a,b+1) - inv(a,b) = (n-1) - 2 p;
  L4  mixed second difference of S = 4n [pi(a) = b] - 4, and the rectangle form
      (= 4n times the discrepancy of the rectangle);
  L5  the closed form of S(a, b);
  L6  duality S_{pi*}(a, b) = -S_pi(a, -b) for pi*(i) = (-1 - pi(i)) mod n;
  L7  parity of inv_3 of every triple is a cut invariant, and
      (n-2) S(a,b) = sum over triples of (3 - 2 inv_3).

Part 2 (C38): sum_a min_b inv(a, b) <= n * floor((n-1)^2/4), exhaustively for
  4 <= n <= AMAX, with the equality set = the n reflections.  (4 <= n <= 12 is
  covered by experiments/torus_cut_scan.c, data/runs/torus_cut_scan/.)

Part 3 (C39): the refuted strategies, each by its witness:
  (a) one cut only: max over lines y of min_b inv(y - b) - floor((n-1)^2/4)
      = 0,0,0,0,1,1,1,3,3,3,5 for n = 3..13 -- the excess is not bounded;
  (b) induction by deleting one or two points fails at n = 8, pi(i) = 5i mod 8;
  (c) averaging over any linear family {(a, c a + d)} fails;
  (d) averaging over the pi-adapted family {(a, pi(a) - t)} fails;
  (e) the separable relaxation min_{a,b} sum_i min(u_i, n-1-v_i) fails.

Usage: python3 checks/check_C37.py [--amax 8] [--idmax 7]
Output: data/runs/check_C37/report.json, report.md.
"""

import argparse
import itertools
import json
import math
import os
import random
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "data", "runs", "check_C37")
VERSION = "check_C37-1.0"


def _coprime(a, b):
    return math.gcd(a, b) == 1


def inv_line(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def line(perm, a, b):
    n = len(perm)
    return [(perm[(a + j) % n] - b) % n for j in range(n)]


def inv_table(perm):
    n = len(perm)
    return [[inv_line(line(perm, a, b)) for b in range(n)] for a in range(n)]


def inv_table_fast(perm):
    """Same table via the C37/L3 increments (used only where L3 is already checked)."""
    n = len(perm)
    pinv = [0] * n
    for i, v in enumerate(perm):
        pinv[v] = i
    tab = [[0] * n for _ in range(n)]
    tab[0][0] = inv_line(list(perm))
    for b in range(n - 1):
        tab[0][b + 1] = tab[0][b] + (n - 1) - 2 * pinv[b]
    for a in range(n - 1):
        pa = perm[a]
        row, nxt = tab[a], tab[a + 1]
        for b in range(n):
            nxt[b] = row[b] + (n - 1) - 2 * ((pa - b) % n)
    return tab


def s_table(perm):
    n = len(perm)
    B = n * (n - 1) // 2
    return [[B - 2 * v for v in row] for row in inv_table(perm)]


def bound(n):
    return ((n - 1) ** 2) // 4


# ---------------------------------------------------------------- part 1

def check_identities(n, sample=None, seed=12345):
    """All C37 lemmas for every permutation of Z_n (or for a sample of them)."""
    B = n * (n - 1) // 2
    if sample is None:
        pool = itertools.permutations(range(n))
    else:
        rnd = random.Random(seed)
        pool = [tuple((h - i) % n for i in range(n)) for h in range(n)]
        pool += [tuple((c * i + d) % n for i in range(n))
                 for c in range(1, n) for d in range(n) if _coprime(c, n)]
        while len(pool) < sample:
            q = list(range(n))
            rnd.shuffle(q)
            pool.append(tuple(q))
    for perm in pool:
        pinv = [0] * n
        for i, v in enumerate(perm):
            pinv[v] = i
        S = s_table(perm)
        iv = [[(B - S[a][b]) // 2 for b in range(n)] for a in range(n)]

        # L1: S = concordant - discordant, threshold
        for a in range(n):
            for b in range(n):
                conc = B - iv[a][b]
                assert S[a][b] == conc - iv[a][b]
        assert (max(max(r) for r in S) >= n // 2) == (min(min(r) for r in iv) <= bound(n))

        # L2: agreements of the two linearizations
        for a in (0, 1 % n):
            for b in (0, 1 % n):
                agree = 0
                for i in range(n):
                    for j in range(i + 1, n):
                        pi_lt = ((i - a) % n) < ((j - a) % n)
                        pv_lt = ((perm[i] - b) % n) < ((perm[j] - b) % n)
                        if pi_lt == pv_lt:
                            agree += 1
                assert agree == B - iv[a][b], (perm, a, b)

        assert inv_table_fast(perm) == iv, perm

        # L3: increments
        for a in range(n):
            for b in range(n):
                w0 = (perm[a] - b) % n
                p = (pinv[b] - a) % n
                assert iv[(a + 1) % n][b] - iv[a][b] == (n - 1) - 2 * w0
                assert iv[a][(b + 1) % n] - iv[a][b] == (n - 1) - 2 * p

        # L4: mixed second difference and the rectangle form
        for a in range(n):
            for b in range(n):
                mix = (S[(a + 1) % n][(b + 1) % n] - S[(a + 1) % n][b]
                       - S[a][(b + 1) % n] + S[a][b])
                assert mix == 4 * n * (1 if perm[a] == b else 0) - 4
        for a1 in range(n + 1):
            for a2 in range(a1, n + 1):
                for b1 in range(n + 1):
                    for b2 in range(b1, n + 1):
                        pts = sum(1 for i in range(a1, a2) if b1 <= perm[i % n] < b2)
                        area = (a2 - a1) * (b2 - b1)
                        lhs = (S[a2 % n][b2 % n] + S[a1 % n][b1 % n]
                               - S[a1 % n][b2 % n] - S[a2 % n][b1 % n])
                        assert lhs * n == 4 * n * (pts * n - area), (perm, a1, a2, b1, b2)

        # L5: closed form
        for a in range(n + 1):
            for b in range(n + 1):
                K = sum(1 for i in range(a) if perm[i % n] < b)
                val = (S[0][0]
                       + 4 * sum(perm[i % n] for i in range(a)) - 2 * a * (n - 1)
                       + 4 * sum(pinv[v % n] for v in range(b)) - 2 * b * (n - 1)
                       + 4 * n * K - 4 * a * b)
                assert val == S[a % n][b % n], (perm, a, b)

        # L6: duality
        star = tuple((-1 - v) % n for v in perm)
        Ss = s_table(star)
        for a in range(n):
            for b in range(n):
                assert Ss[a][b] == -S[a][(-b) % n]

        # L7: triples
        for T in itertools.combinations(range(n), 3):
            par = None
            for a in range(n):
                for b in range(n):
                    w = line(perm, a, b)
                    pos = sorted(((i - a) % n, (perm[i] - b) % n) for i in T)
                    v = [y for _, y in pos]
                    i3 = sum(1 for x in range(3) for y in range(x + 1, 3) if v[x] > v[y])
                    if par is None:
                        par = i3 % 2
                    assert i3 % 2 == par, (perm, T, a, b)
        for a in range(n):
            for b in range(n):
                tot = 0
                for T in itertools.combinations(range(n), 3):
                    pos = sorted(((i - a) % n, (perm[i] - b) % n) for i in T)
                    v = [y for _, y in pos]
                    i3 = sum(1 for x in range(3) for y in range(x + 1, 3) if v[x] > v[y])
                    tot += 3 - 2 * i3
                assert tot == (n - 2) * S[a][b], (perm, a, b)
    return True


def check_reflection_value(nmax):
    """I(pi_h) = floor((n-1)^2/4) for every reflection, 4 <= n <= nmax."""
    for n in range(4, nmax + 1):
        for h in range(n):
            perm = tuple((h - i) % n for i in range(n))
            best = min(min(row) for row in inv_table_fast(perm))
            assert best == bound(n), (n, h, best)
    return True


# ---------------------------------------------------------------- part 2

def check_row_sum(n):
    """C38: sum_a min_b inv(a, b) <= n * floor((n-1)^2/4); equality set."""
    target = n * bound(n)
    worst, eq = -1, []
    for perm in itertools.permutations(range(n)):
        iv = inv_table_fast(perm)
        tot = sum(min(row) for row in iv)
        if tot > worst:
            worst, eq = tot, [perm]
        elif tot == worst:
            eq.append(perm)
    reflections = {tuple((h - i) % n for i in range(n)) for h in range(n)}
    return {"n": n, "max": worst, "target": target, "ok": worst <= target,
            "equality_count": len(eq), "equality_is_reflections": set(eq) == reflections}


# ---------------------------------------------------------------- part 3

def min_shift_inv(y):
    """min over value shifts b of inv(y - b mod n)."""
    n = len(y)
    pos = [0] * n
    for i, v in enumerate(y):
        pos[v] = i
    cur = inv_line(list(y))
    best = cur
    for b in range(n - 1):
        cur += (n - 1) - 2 * pos[b]
        best = min(best, cur)
    return best


def contract(perm, dels):
    n = len(perm)
    vals = {perm[i] for i in dels}
    keep = [i for i in range(n) if i not in dels]
    vs = sorted(v for v in range(n) if v not in vals)
    rank = {v: k for k, v in enumerate(vs)}
    return tuple(rank[perm[i]] for i in keep)


def I_of(perm):
    return min(min(row) for row in inv_table_fast(perm))


def check_negative(idmax):
    res = {}

    # (a) one cut only: exhaustive for 3..idmax, witnesses for 10..13
    exc = {}
    for n in range(3, idmax + 1):
        best = max(min_shift_inv(y) for y in itertools.permutations(range(n)))
        exc[n] = best - bound(n)
    wit = {10: (0, 8, 6, 5, 4, 3, 2, 1, 9, 7),
           11: (0, 9, 7, 6, 3, 5, 4, 2, 1, 10, 8),
           12: (0, 2, 10, 8, 7, 6, 5, 4, 3, 11, 9, 1),
           13: (0, 9, 11, 7, 6, 5, 4, 3, 2, 1, 12, 10, 8)}
    wit_exc = {n: min_shift_inv(y) - bound(n) for n, y in wit.items()}
    res["one_cut_excess_exhaustive"] = exc
    res["one_cut_excess_witnesses"] = wit_exc
    res["one_cut_unbounded"] = max(wit_exc.values()) >= 5

    # (b) deletion induction, n = 8, pi(i) = 5 i mod 8
    n = 8
    perm = tuple((5 * i) % n for i in range(n))
    I = I_of(perm)
    d1 = [I_of(contract(perm, {i})) for i in range(n)]
    d2 = [I_of(contract(perm, {i, j})) for i in range(n) for j in range(i + 1, n)]
    res["deletion"] = {
        "perm": list(perm), "I": I,
        "max_I_after_1_deletion": max(d1),
        "step_needed_1": bound(n) - bound(n - 1),
        "ok_1": I - max(d1) <= bound(n) - bound(n - 1),
        "max_I_after_2_deletions": max(d2),
        "step_needed_2": bound(n) - bound(n - 2),
        "ok_2": I - max(d2) <= bound(n) - bound(n - 2)}

    # (c) linear families {(a, c a + d)}
    lin = []
    for n, perm in ((5, (0, 2, 1, 4, 3)), (7, (0, 2, 4, 1, 6, 5, 3))):
        S = s_table(perm)
        best = max(sum(S[a][(c * a + d) % n] for a in range(n))
                   for c in range(n) for d in range(n))
        lin.append({"n": n, "perm": list(perm), "best_sum": best,
                    "target": n * (n // 2), "fails": best < n * (n // 2)})
    res["linear_family"] = lin

    # (d) pi-adapted family {(a, pi(a) - t)}
    adp = []
    for n, perm in ((5, (0, 1, 4, 3, 2)), (7, (0, 2, 1, 6, 5, 4, 3))):
        S = s_table(perm)
        best = max(sum(S[a][(perm[a] - t) % n] for a in range(n)) for t in range(n))
        adp.append({"n": n, "perm": list(perm), "best_sum": best,
                    "target": n * (n // 2), "fails": best < n * (n // 2)})
    res["adapted_family"] = adp

    # (e) separable relaxation: exhaustive excess for 4 <= n <= 7, witness at n = 8
    def minG(perm):
        n = len(perm)
        return min(sum(min((i - a) % n, n - 1 - ((perm[i] - b) % n)) for i in range(n))
                   for a in range(n) for b in range(n))
    sep = {n: max(minG(p) for p in itertools.permutations(range(n))) - bound(n)
           for n in range(4, 8)}
    w8 = tuple((5 * i) % 8 for i in range(8))
    res["separable_excess_4_7"] = sep
    res["separable_witness_n8"] = {"perm": list(w8), "min_G": minG(w8), "bound": bound(8)}
    res["separable_fails"] = (max(sep.values()) > 0
                              and minG(w8) - bound(8) == 4)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amax", type=int, default=8, help="exhaustive C38 up to this n")
    ap.add_argument("--idmax", type=int, default=6, help="exhaustive identities up to this n")
    ap.add_argument("--idsample_max", type=int, default=8, help="sampled identities up to this n")
    ap.add_argument("--idsample", type=int, default=300, help="sample size per n")
    ap.add_argument("--onecut_max", type=int, default=9,
                    help="exhaustive one-cut excess up to this n")
    args = ap.parse_args()

    t0 = time.time()
    rep = {"version": VERSION, "amax": args.amax, "idmax": args.idmax}

    ident, ident_sampled = {}, {}
    for n in range(4, args.idmax + 1):
        t = time.time()
        ident[n] = check_identities(n)
        print("identities n=%d exhaustive: PASS (%.1f s)" % (n, time.time() - t))
    for n in range(args.idmax + 1, args.idsample_max + 1):
        t = time.time()
        ident_sampled[n] = check_identities(n, sample=args.idsample)
        print("identities n=%d sampled(%d): PASS (%.1f s)"
              % (n, args.idsample, time.time() - t))
    rep["identities_exhaustive"] = ident
    rep["identities_sampled"] = ident_sampled
    rep["identities_sample_size"] = args.idsample
    rep["reflection_value_4_20"] = check_reflection_value(20)
    print("reflections I(pi_h) = floor((n-1)^2/4), 4 <= n <= 20: PASS")

    rows = []
    for n in range(4, args.amax + 1):
        r = check_row_sum(n)
        rows.append(r)
        print("C38 n=%d: max=%d target=%d ok=%s eq=%d reflections=%s"
              % (n, r["max"], r["target"], r["ok"], r["equality_count"],
                 r["equality_is_reflections"]))
    rep["row_sum"] = rows

    rep["negative"] = check_negative(args.onecut_max)
    print("negative results:", json.dumps(rep["negative"], sort_keys=True)[:400], "...")

    ok = (all(ident.values()) and all(ident_sampled.values())
          and rep["reflection_value_4_20"]
          and all(r["ok"] and r["equality_is_reflections"] for r in rows)
          and rep["negative"]["one_cut_unbounded"]
          and not rep["negative"]["deletion"]["ok_1"]
          and not rep["negative"]["deletion"]["ok_2"]
          and all(r["fails"] for r in rep["negative"]["linear_family"])
          and all(r["fails"] for r in rep["negative"]["adapted_family"])
          and rep["negative"]["separable_fails"])
    rep["seconds"] = round(time.time() - t0, 1)
    rep["result"] = "PASS" if ok else "FAIL"

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(rep, f, indent=1, sort_keys=True)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# check_C37 — %s\n\n" % rep["result"])
        f.write("Версия %s, %.1f с.\n\n" % (VERSION, rep["seconds"]))
        f.write("- тождества C37 (леммы 1–7): все pi при 4 <= n <= %d (исчерпывающе); "
                "выборка %d перестановок (отражения, аффинные, случайные, seed 12345) "
                "при %d <= n <= %d\n"
                % (args.idmax, args.idsample, args.idmax + 1, args.idsample_max))
        f.write("- `I(pi_h) = floor((n-1)^2/4)` на всех отражениях при 4 <= n <= 20\n")
        f.write("- C38 `sum_a min_b inv <= n floor((n-1)^2/4)`: все pi при 4 <= n <= %d, "
                "равенство ровно на n отражениях\n" % args.amax)
        f.write("- C39 (отрицательные результаты): избыток одного разреза %s; "
                "индукция удалением 1 и 2 точек не проходит; усреднения по линейным "
                "и по pi-адаптированным семействам не проходят; сепарабельная "
                "релаксация не проходит\n"
                % json.dumps(rep["negative"]["one_cut_excess_exhaustive"]))
    print(rep["result"], "%.1f s" % rep["seconds"])
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
