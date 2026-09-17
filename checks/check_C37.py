"""Checker for C37/C38 (session 9): H13-I as a switching problem.

Independent implementation (nothing is imported from experiments/ or
constructions/; only the certified distance tables are *not* needed here, the
statement is purely combinatorial -- inversions of lines of a toric class).

Parts:
  L1  Lemma 1: inv <= floor((n-1)^2/4)  <=>  N - 2 inv >= floor(n/2),
      N = n(n-1)/2, for 2 <= n <= 200 and all inv.
  H   H13-I itself: max over the n^2 double cuts of S = N - 2 inv is at least
      floor(n/2) -- exhaustively over all toric class representatives
      (pi(0) = 0) for 4 <= n <= HMAX, and on families (reflections, affine,
      random, block words) for n up to 40.  Also re-derives max_pi I(pi) =
      floor((n-1)^2/4) with equality exactly on the n reflections (C33).
  L2  Lemma 2: chi(q,s) = chi(0,0) switched by A_q xor B_s, and one step of a
      cut is a one-vertex switching -- all pi, all cuts, 4 <= n <= 7.
  T   Theorem 4: every +-1 signing of K_n has a switching with S >= floor(n/2),
      and the all-minus class attains floor(n/2) exactly.  All 2^C(n,2)
      signings for n = 4, 5; signings of permutations for n = 6, 7, 8; random
      signings for 6 <= n <= 9.
  L7  Lemma 7 (identity S = (4/n)K + (2/n^2)T - 1 - n/2), integer form
      2 n^2 S = 2 n K2 + T2 - 2 n^2 - n^3 with K2 = 4K, T2 = 4T.
  P6  Proposition 6 / C38: permutations with no stable cut (every cut has an
      element in more than floor((n-1)/2) inversions).  Exhaustive for
      4 <= n <= 8: none for n <= 7, exactly pi(i) = 3i and pi(i) = 5i for n = 8.
  R   Corollary 5: for reflections every cut gives the line
      (c, c-1, ..., 0, n-1, ..., c+1), inv = C(c+1,2) + C(n-1-c,2), the minimum
      over cuts is exactly floor((n-1)^2/4) and the all-minus signing occurs;
      4 <= n <= 24, all h.
  V   The nu-average family (shifts of the cut along the permutation's own
      points) does not certify H13-I: for n = 5 and n = 7 the maximum over
      (a,b) of the average of S over the n cuts (j+a, pi(j)+b) is < floor(n/2).

Usage: python3 checks/check_C37.py [--hmax 8]
Output: data/runs/check_C37/report.json, report.md
"""

import argparse
import itertools
import json
import os
import random
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
VERSION = "check_C37-1.0"
OUT = os.path.join(ROOT, "data", "runs", "check_C37")


# ---------------------------------------------------------------- primitives

def line(pi, q, s):
    n = len(pi)
    return [(pi[(q + j) % n] - s) % n for j in range(n)]


def inv_and_degrees(w):
    """inversions of w and, per line index, the number of inversions it is in"""
    n = len(w)
    r = [0] * n
    inv = 0
    for i in range(n):
        wi = w[i]
        for j in range(i + 1, n):
            if wi > w[j]:
                r[i] += 1
                r[j] += 1
                inv += 1
    return inv, r


def reps(n):
    """toric class representatives: pi(0) = 0 (value rotation fixes it)"""
    for tail in itertools.permutations(range(1, n)):
        yield (0,) + tail


def reflections(n):
    return [tuple((h - i) % n for i in range(n)) for h in range(n)]


# ------------------------------------------------------------------ part L1

def part_L1(log):
    ok = True
    for n in range(2, 201):
        N = n * (n - 1) // 2
        bnd = (n - 1) ** 2 // 4
        for inv in range(0, N + 1):
            if (inv <= bnd) != (N - 2 * inv >= n // 2):
                ok = False
        if (N - 2 * bnd) != n // 2:
            ok = False
    log("L1: equivalence inv <= floor((n-1)^2/4) <=> S >= floor(n/2), 2<=n<=200: %s"
        % ("PASS" if ok else "FAIL"))
    return ok, {"range": "2<=n<=200"}


# ------------------------------------------------------------------- part H

def best_cut_S(pi):
    n = len(pi)
    N = n * (n - 1) // 2
    best = None
    besti = None
    for q in range(n):
        for s in range(n):
            inv, _ = inv_and_degrees(line(pi, q, s))
            if besti is None or inv < besti:
                besti = inv
                best = N - 2 * inv
    return best, besti


def part_H(hmax, log):
    ok = True
    rows = []
    for n in range(4, hmax + 1):
        thr = n // 2
        bnd = (n - 1) ** 2 // 4
        maxI = -1
        argmaxI = []
        cnt = 0
        t0 = time.time()
        for pi in reps(n):
            cnt += 1
            S, I = best_cut_S(pi)
            if S < thr or I > bnd:
                ok = False
                log("  H FAIL n=%d pi=%s S=%d I=%d" % (n, pi, S, I))
            if I > maxI:
                maxI, argmaxI = I, [pi]
            elif I == maxI:
                argmaxI.append(pi)
        refl = set()
        for r in reflections(n):
            # representative of the reflection's toric class with pi(0)=0
            refl.add(tuple((v - r[0]) % n for v in r))
        eq_refl = set(argmaxI) == refl
        rows.append({"n": n, "reps": cnt, "max_I": maxI, "bound": bnd,
                     "argmax_are_reflections": eq_refl,
                     "n_argmax": len(argmaxI), "sec": round(time.time() - t0, 1)})
        log("H: n=%d reps=%d max_I=%d (bound %d) argmax=%d classes, "
            "exactly the reflections: %s [%.1fs]"
            % (n, cnt, maxI, bnd, len(argmaxI), eq_refl, time.time() - t0))
        if maxI != bnd or not eq_refl:
            ok = False
    # families for larger n
    rnd = random.Random(20260917)
    fam_rows = []
    for n in list(range(4, 21)) + [24, 30, 33, 40]:
        thr = n // 2
        bnd = (n - 1) ** 2 // 4
        fam = []
        fam += [("refl%d" % h, tuple((h - i) % n for i in range(n))) for h in range(n)]
        for a in range(1, n):
            if all(a * k % n != 0 for k in range(1, n)) and \
               len({a * i % n for i in range(n)}) == n:
                fam.append(("aff%d" % a, tuple((a * i) % n for i in range(n))))
        for k in range(5):
            p = list(range(n))
            rnd.shuffle(p)
            fam.append(("rand%d" % k, tuple(p)))
        b = n // 2
        fam.append(("block", tuple(list(range(b, n)) + list(range(b)))))
        worst = None
        for name, pi in fam:
            S, I = best_cut_S(pi)
            if S < thr or I > bnd:
                ok = False
                log("  H FAIL (family) n=%d %s S=%d I=%d" % (n, name, S, I))
            if worst is None or S < worst[0]:
                worst = (S, name, I)
        fam_rows.append({"n": n, "family_size": len(fam), "min_S": worst[0],
                         "at": worst[1], "I": worst[2], "thr": thr, "bound": bnd})
    log("H: families (reflections, affine, random, block) n<=40: min S over family "
        "always >= floor(n/2): %s" % ("PASS" if ok else "FAIL"))
    return ok, {"exhaustive": rows, "families": fam_rows}


# ------------------------------------------------------------------ part L2

def signing(pi, q, s):
    """chi[u][v] for the double cut (q,s); vertices are circle positions u."""
    n = len(pi)
    place = {}
    for j in range(n):
        u = (q + j) % n
        place[u] = (j, (pi[u] - s) % n)
    chi = [[0] * n for _ in range(n)]
    for u in range(n):
        for v in range(n):
            if u == v:
                continue
            ju, yu = place[u]
            jv, yv = place[v]
            chi[u][v] = 1 if (ju - jv) * (yu - yv) > 0 else -1
    return chi


def switched(chi, U, n):
    c = [row[:] for row in chi]
    for u in range(n):
        for v in range(n):
            if u != v and ((u in U) != (v in U)):
                c[u][v] = -c[u][v]
    return c


def part_L2(log):
    ok = True
    for n in range(4, 8):
        for pi in reps(n):
            base = signing(pi, 0, 0)
            for q in range(n):
                for s in range(n):
                    A = set(range(q))
                    B = {i for i in range(n) if pi[i] < s}
                    if switched(base, A ^ B, n) != signing(pi, q, s):
                        ok = False
                    # one-step: cut (q,s) -> (q+1,s) switches vertex q
                    if switched(signing(pi, q, s), {q}, n) != \
                            signing(pi, (q + 1) % n, s):
                        ok = False
                    iv = pi.index(s % n)
                    if switched(signing(pi, q, s), {iv}, n) != \
                            signing(pi, q, (s + 1) % n):
                        ok = False
    log("L2: cut family = switchings by A_q xor B_s, one step = one vertex, "
        "all pi and cuts 4<=n<=7: %s" % ("PASS" if ok else "FAIL"))
    return ok, {"range": "4<=n<=7"}


# ------------------------------------------------------------------- part T

def S_of(chi, n):
    return sum(chi[u][v] for u in range(n) for v in range(u + 1, n))


def best_switch_S(chi, n):
    best = None
    for mask in range(1 << (n - 1)):
        U = {i for i in range(n - 1) if mask >> i & 1}
        S = S_of(switched(chi, U, n), n)
        if best is None or S > best:
            best = S
    return best


def part_T(log):
    ok = True
    detail = {}
    # (a) all signings for n = 4, 5
    for n in (4, 5):
        pairs = [(u, v) for u in range(n) for v in range(u + 1, n)]
        thr = n // 2
        worst = None
        for mask in range(1 << len(pairs)):
            chi = [[0] * n for _ in range(n)]
            for b, (u, v) in enumerate(pairs):
                x = 1 if mask >> b & 1 else -1
                chi[u][v] = chi[v][u] = x
            b = best_switch_S(chi, n)
            if b < thr:
                ok = False
            if worst is None or b < worst:
                worst = b
        detail["all_signings_n%d" % n] = {"min_max_S": worst, "thr": thr,
                                          "signings": 1 << len(pairs)}
        log("T: n=%d, all %d signings: min over signings of max over switchings "
            "S = %d (need >= %d)" % (n, 1 << len(pairs), worst, thr))
        if worst != thr:
            ok = False
    # (b) signings coming from permutations, n = 6,7,8
    for n in (6, 7, 8):
        thr = n // 2
        worst = None
        gaps = {}
        for pi in reps(n):
            chi = signing(pi, 0, 0)
            b = best_switch_S(chi, n)
            c, _ = best_cut_S(pi)
            if b < thr or c < thr:
                ok = False
            g = b - c
            gaps[g] = gaps.get(g, 0) + 1
            if worst is None or b < worst:
                worst = b
        detail["perm_signings_n%d" % n] = {"min_max_switch_S": worst,
                                           "thr": thr, "gap_hist": gaps}
        log("T: n=%d, permutation signings: min max_switch S = %d (need >= %d); "
            "hist of (max_switch - max_cut) = %s" % (n, worst, thr, gaps))
    # (c) random signings
    rnd = random.Random(4242)
    for n in range(6, 10):
        thr = n // 2
        worst = None
        for _ in range(300):
            chi = [[0] * n for _ in range(n)]
            for u in range(n):
                for v in range(u + 1, n):
                    x = rnd.choice((-1, 1))
                    chi[u][v] = chi[v][u] = x
            b = best_switch_S(chi, n)
            if b < thr:
                ok = False
            if worst is None or b < worst:
                worst = b
        detail["random_signings_n%d" % n] = {"trials": 300, "min_max_S": worst,
                                             "thr": thr}
    # (d) tightness of the all-minus class
    for n in range(2, 41):
        N = n * (n - 1) // 2
        best = max(2 * u * (n - u) - N for u in range(n + 1))
        if best != n // 2:
            ok = False
    log("T: all-minus class maximum = floor(n/2) exactly, 2<=n<=40; theorem: %s"
        % ("PASS" if ok else "FAIL"))
    return ok, detail


# ------------------------------------------------------------------ part L7

def part_L7(hmax, log):
    """2 n^2 S = 2 n K2 + T2 - 2 n^2 - n^3, K2 = 4K, T2 = 4T (integers)."""
    ok = True
    for n in range(4, min(hmax, 8) + 1):
        N = n * (n - 1) // 2

        def s2(t):
            return 2 * (t % n) - n

        for pi in reps(n):
            T2 = sum(s2(j - i) * s2(pi[j] - pi[i])
                     for i in range(n) for j in range(n))
            for q in range(n):
                for s in range(n):
                    K2 = sum(s2(i - q) * s2(pi[i] - s) for i in range(n))
                    inv, _ = inv_and_degrees(line(pi, q, s))
                    S = N - 2 * inv
                    if 2 * n * n * S != 2 * n * K2 + T2 - 2 * n * n - n ** 3:
                        ok = False
    log("L7: identity S = (4/n)K + (2/n^2)T - 1 - n/2, all pi and cuts "
        "4<=n<=%d: %s" % (min(hmax, 8), "PASS" if ok else "FAIL"))
    # mean over cuts for reflections
    rows = []
    for n in range(4, 13):
        pi = tuple((-i) % n for i in range(n))
        tot = 0
        for q in range(n):
            for s in range(n):
                inv, _ = inv_and_degrees(line(pi, q, s))
                tot += inv
        mean_noninv = n * (n - 1) / 2 - tot / n ** 2
        pred = (n * n - 1) / 6
        rows.append({"n": n, "mean_noninv": mean_noninv, "(n^2-1)/6": pred,
                     "floor(n^2/4)": n * n // 4})
        if abs(mean_noninv - pred) > 1e-9:
            ok = False
    log("L7: mean number of non-inversions over all n^2 cuts of a reflection "
        "= (n^2-1)/6 < floor(n^2/4), 4<=n<=12: PASS")
    return ok, {"reflection_mean": rows}


# ------------------------------------------------------------------ part P6

def part_P6(log):
    ok = True
    detail = []
    for n in range(4, 9):
        thr = (n - 1) // 2
        bad = []
        for pi in reps(n):
            stable = False
            for q in range(n):
                for s in range(n):
                    _, r = inv_and_degrees(line(pi, q, s))
                    if max(r) <= thr:
                        stable = True
                        break
                if stable:
                    break
            if not stable:
                bad.append(pi)
        detail.append({"n": n, "no_stable_cut": [list(p) for p in bad],
                       "count": len(bad)})
        log("P6: n=%d, representatives with no stable cut: %d %s"
            % (n, len(bad), [list(p) for p in bad] if bad else ""))
        if n <= 7 and bad:
            ok = False
        if n == 8:
            want = {tuple((3 * i) % 8 for i in range(8)),
                    tuple((5 * i) % 8 for i in range(8))}
            if set(bad) != want:
                ok = False
    log("P6: no stable cut for 4<=n<=7: none; n=8: exactly pi(i)=3i and 5i: %s"
        % ("PASS" if ok else "FAIL"))
    return ok, detail


# ------------------------------------------------------------------- part V

def part_V(log):
    """the nu-family (cuts (j+a, pi(j)+b) averaged over j) is not enough"""
    ok = True
    rows = []
    for n in (5, 7):
        thr = n // 2
        worst = None
        N = n * (n - 1) // 2
        for pi in reps(n):
            St = [[N - 2 * inv_and_degrees(line(pi, q, s))[0] for s in range(n)]
                  for q in range(n)]
            best = None
            for a in range(n):
                for b in range(n):
                    tot = sum(St[(j + a) % n][(pi[j] + b) % n] for j in range(n))
                    if best is None or tot > best:
                        best = tot
            if worst is None or best < worst[0]:
                worst = (best, pi)
        rows.append({"n": n, "min_over_pi of max_ab n*Phi": worst[0],
                     "n*thr": n * thr, "at": list(worst[1])})
        log("V: n=%d: min over pi of max over (a,b) of sum_j S(j+a, pi(j)+b) = %d"
            " < n*floor(n/2) = %d at pi=%s"
            % (n, worst[0], n * thr, list(worst[1])))
        if worst[0] >= n * thr:
            ok = False
    return ok, rows


# ------------------------------------------------------------------- part R

def part_R(log):
    """Corollary 5: for a reflection every cut gives the line
    (c, c-1, ..., 0, n-1, ..., c+1) with inv = C(c+1,2) + C(n-1-c,2), the
    minimum over cuts is exactly floor((n-1)^2/4), and c = n-1 gives the
    all-minus signing."""
    from math import comb
    ok = True
    rows = []
    for n in range(4, 25):
        for h in range(n):
            pi = [(h - i) % n for i in range(n)]
            vals = []
            for q in range(n):
                for s in range(n):
                    w = [(pi[(q + j) % n] - s) % n for j in range(n)]
                    c = w[0]
                    if w != [(c - j) % n for j in range(n)]:
                        ok = False
                    inv, _ = inv_and_degrees(w)
                    if inv != comb(c + 1, 2) + comb(n - 1 - c, 2):
                        ok = False
                    vals.append((inv, c))
            mn = min(v[0] for v in vals)
            full = max(v[0] for v in vals)
            if mn != (n - 1) ** 2 // 4 or full != n * (n - 1) // 2:
                ok = False
        rows.append({"n": n, "min_inv": mn, "bound": (n - 1) ** 2 // 4})
    log("R: reflections -- line form, inv = C(c+1,2)+C(n-1-c,2), min over cuts "
        "= floor((n-1)^2/4), all-minus signing present, 4<=n<=24, all h: %s"
        % ("PASS" if ok else "FAIL"))
    return ok, rows


# --------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hmax", type=int, default=8)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    lines = []

    def log(msg):
        print(msg)
        lines.append(msg)

    t0 = time.time()
    res = {}
    allok = True
    for name, fn in (("L1", lambda: part_L1(log)),
                     ("L2", lambda: part_L2(log)),
                     ("T", lambda: part_T(log)),
                     ("L7", lambda: part_L7(args.hmax, log)),
                     ("P6", lambda: part_P6(log)),
                     ("V", lambda: part_V(log)),
                     ("R", lambda: part_R(log)),
                     ("H", lambda: part_H(args.hmax, log))):
        ok, d = fn()
        res[name] = {"ok": ok, "detail": d}
        allok = allok and ok
    res["version"] = VERSION
    res["hmax"] = args.hmax
    res["seconds"] = round(time.time() - t0, 1)
    res["PASS"] = allok
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(res, f, indent=1, sort_keys=True, default=str)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# check_C37 — отчёт\n\n%s, hmax=%d, %.1f с.\n\n```text\n%s\n```\n"
                % (VERSION, args.hmax, res["seconds"], "\n".join(lines)))
    print("PASS" if allok else "FAIL", "%.1f s" % res["seconds"])
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
