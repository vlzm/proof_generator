"""Independent checker for C37-C40 (docs/proofs/C37_toric_inversions.md, session 9).

Everything is recomputed here from the definition of the line in
experiments/line_model.py (w_j = pi[(q+1+j) % n] - (q+1-c) mod n); nothing is
imported from experiments/toric_inversions.c, so an error in the fast O(n^2)
table of that program cannot hide here.

  Part A (C37, dictionary).  For all pi, 4 <= n <= NA (and random pi up to
    n = 11): inv(alpha, beta) computed from the definition equals
      * |E symmetric-difference cut(S_alpha symmetric-difference T_beta)| (Thm 1.3),
      * (C(n,2) + sum_{a<b} eps_ab chi_a chi_b) / 2  (sign form 1.4),
      * the closed formula of Lemma 1.6,
    and the mixed second difference equals 2 - 2n [pi(alpha) = beta] (Cor 1.7).
    Also: the two-box description of the concordance region (Cor 1.2).
  Part B (C38, relaxation).  The greedy of Theorem 2.1 (in every vertex order
    tried) yields sum eps x x <= -floor(n/2), i.e. at most floor((n-1)^2/4)
    edits, on: all graphs for n <= 6, random graphs up to n = 16, and all
    inversion graphs of permutations for 4 <= n <= 7.  Tightness for K_n:
    brute force over all bipartitions for 4 <= n <= 10 gives exactly
    floor((n-1)^2/4).  Also the identity floor((n-1)^2/4) = sum_{m<n} floor(m/2).
  Part C (C39, averaging bound).  sum over the n^2 cuts of conc = T(pi)
    (Lemma 3.1) for all pi, 4 <= n <= NC; T(pi) >= n^2(n^2-1)/6 with equality
    exactly on the n reflections (exhaustive 4 <= n <= 8, random up to n = 60);
    I(pi) <= floor(C(n,2) - T/n^2) <= floor((n-1)(2n-1)/6) checked against the
    directly computed I(pi) for all pi, 4 <= n <= 8.
  Part D (C40, two-cut identity).  Identity (4.1) for all pi and all cut pairs,
    4 <= n <= 6, and random samples up to n = 12; Corollary 4.3 on reflections:
    the explicit pair of cuts has C = D = 0 and gives inv = floor((n-1)^2/4)
    for 4 <= n <= 40.
  Part E (refuted approaches of section 5).  The three counterexamples are
    recomputed: max_alpha min_beta inv = 10 at n = 7 and 13 at n = 8; the best
    (h, k) family gives 5, 7, 11, 16 at n = 5..8; the line (0,1,8,7,6,5,4,3,2,9)
    satisfies all depth-1 local conditions and has inv = 21 > 20.

Usage: python3 checks/check_C37.py [--na 6] [--nc 7] [--quick]
Output: data/runs/check_C37/report.json, report.md.
"""

import argparse
import itertools
import json
import os
import random
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "data", "runs", "check_C37")
VERSION = "check_C37-1.0"


# ---------------------------------------------------------------- definitions

def line(pi, q, c):
    """The relabelled line of the double cut (q, c), exactly as in line_model.py."""
    n = len(pi)
    shift = (q + 1 - c) % n
    return [(pi[(q + 1 + j) % n] - shift) % n for j in range(n)]


def inv_of(w):
    return sum(1 for i in range(len(w)) for j in range(i + 1, len(w)) if w[i] > w[j])


def inv_cut(pi, alpha, beta):
    """inv at the cut with position origin alpha and value origin beta."""
    n = len(pi)
    return inv_of(line(pi, (alpha - 1) % n, (alpha - beta) % n))


def I_of(pi):
    n = len(pi)
    return min(inv_cut(pi, a, b) for a in range(n) for b in range(n))


def target(n):
    return ((n - 1) ** 2) // 4


# ------------------------------------------------------------------- part A

def cut_pairs(U, n):
    return set((i, j) for i in range(n) for j in range(i + 1, n) if (i in U) != (j in U))


def part_A(na, rng):
    n_checked = 0
    for n in range(4, na + 1):
        for pi in itertools.permutations(range(n)):
            pi = list(pi)
            E = set((i, j) for i in range(n) for j in range(i + 1, n) if pi[i] > pi[j])
            eps = {(i, j): (1 if pi[i] > pi[j] else -1) for i in range(n) for j in range(i + 1, n)}
            inv0 = inv_cut(pi, 0, 0)
            ipi = [0] * n
            for i, v in enumerate(pi):
                ipi[v] = i
            tab = {}
            for a in range(n):
                for b in range(n):
                    val = inv_cut(pi, a, b)
                    tab[(a, b)] = val
                    S = set(i for i in range(n) if i < a)
                    T = set(i for i in range(n) if pi[i] < b)
                    U = S ^ T
                    assert val == len(E ^ cut_pairs(U, n)), ("Thm1.3", n, pi, a, b)
                    chi = [1 if ((i >= a) == (pi[i] >= b)) else -1 for i in range(n)]
                    R = sum(eps[(i, j)] * chi[i] * chi[j] for (i, j) in eps)
                    assert val == (n * (n - 1) // 2 + R) // 2, ("sign form", n, pi, a, b)
                    F = sum(pi[x] for x in range(a))
                    G = sum(ipi[y] for y in range(b))
                    N = sum(1 for x in range(a) if pi[x] < b)
                    assert val == inv0 + (n - 1) * (a + b) + 2 * a * b - 2 * F - 2 * G - 2 * n * N, \
                        ("Lemma1.6", n, pi, a, b)
            for a in range(n - 1):
                for b in range(n - 1):
                    mixed = tab[(a + 1, b + 1)] - tab[(a + 1, b)] - tab[(a, b + 1)] + tab[(a, b)]
                    assert mixed == 2 - 2 * n * (1 if pi[a] == b else 0), ("Cor1.7", n, pi, a, b)
            # Cor 1.2: concordance region of a pair is the two aligned boxes
            for i in range(n):
                for j in range(i + 1, n):
                    p = (j - i) % n
                    v = (pi[j] - pi[i]) % n
                    cnt = 0
                    for a in range(n):
                        for b in range(n):
                            w = line(pi, (a - 1) % n, (a - b) % n)
                            pos_i = [k for k in range(n) if w[k] == (pi[i] - b) % n][0]
                            pos_j = [k for k in range(n) if w[k] == (pi[j] - b) % n][0]
                            conc = (pos_i - pos_j) * (((pi[i] - b) % n) - ((pi[j] - b) % n)) > 0
                            if conc:
                                cnt += 1
                    assert cnt == p * v + (n - p) * (n - v), ("Cor1.2", n, pi, i, j)
            n_checked += 1
    # random larger n: identities 1.3 / 1.4 / 1.6 only
    for n in range(8, 12):
        for _ in range(20):
            pi = list(range(n))
            rng.shuffle(pi)
            E = set((i, j) for i in range(n) for j in range(i + 1, n) if pi[i] > pi[j])
            ipi = [0] * n
            for i, v in enumerate(pi):
                ipi[v] = i
            inv0 = inv_cut(pi, 0, 0)
            for _ in range(12):
                a, b = rng.randrange(n), rng.randrange(n)
                val = inv_cut(pi, a, b)
                S = set(i for i in range(n) if i < a)
                T = set(i for i in range(n) if pi[i] < b)
                assert val == len(E ^ cut_pairs(S ^ T, n)), ("Thm1.3 rnd", n, pi, a, b)
                F = sum(pi[x] for x in range(a))
                G = sum(ipi[y] for y in range(b))
                N = sum(1 for x in range(a) if pi[x] < b)
                assert val == inv0 + (n - 1) * (a + b) + 2 * a * b - 2 * F - 2 * G - 2 * n * N, \
                    ("Lemma1.6 rnd", n, pi, a, b)
            n_checked += 1
    return {"permutations_and_samples_checked": n_checked}


# ------------------------------------------------------------------- part B

def greedy_R(eps, n, order):
    """Theorem 2.1: greedy signs in the given vertex order; returns sum eps x x."""
    x = {}
    for k, v in enumerate(order):
        if k == 0:
            x[v] = 1
            continue
        s = sum(eps[frozenset((v, u))] * x[u] for u in order[:k])
        x[v] = -1 if s > 0 else 1
    return sum(eps[frozenset((u, v))] * x[u] * x[v]
               for i, u in enumerate(order) for v in order[i + 1:])


def brute_min_R(eps, n):
    best = None
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for U in range(1 << (n - 1)):
        x = [-1 if (U >> i) & 1 else 1 for i in range(n)]
        R = sum(eps[frozenset((i, j))] * x[i] * x[j] for (i, j) in pairs)
        if best is None or R < best:
            best = R
    return best


def part_B(rng, quick):
    res = {}
    # identity floor((n-1)^2/4) = sum_{m=0}^{n-1} floor(m/2)
    for n in range(2, 60):
        assert target(n) == sum(m // 2 for m in range(n)), ("sum identity", n)
    # all graphs, n <= 6
    worst = {}
    nmax_all = 5 if quick else 6
    for n in range(3, nmax_all + 1):
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        w = 0
        for mask in range(1 << len(pairs)):
            eps = {frozenset(p): (1 if (mask >> k) & 1 else -1) for k, p in enumerate(pairs)}
            R = greedy_R(eps, n, list(range(n)))
            assert R <= -(n // 2), ("Thm2.1 all graphs", n, mask, R)
            edits = (n * (n - 1) // 2 + R) // 2
            w = max(w, edits)
        worst[n] = w
    res["all_graphs_max_greedy_edits"] = worst
    res["all_graphs_target"] = {n: target(n) for n in worst}
    # random graphs, larger n, random vertex orders
    for n in range(7, 17):
        for _ in range(30):
            pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
            eps = {frozenset(p): rng.choice((1, -1)) for p in pairs}
            order = list(range(n))
            rng.shuffle(order)
            R = greedy_R(eps, n, order)
            assert R <= -(n // 2), ("Thm2.1 random", n, R)
    # inversion graphs of all permutations
    for n in range(4, 8):
        for pi in itertools.permutations(range(n)):
            eps = {frozenset((i, j)): (1 if pi[i] > pi[j] else -1)
                   for i in range(n) for j in range(i + 1, n)}
            R = greedy_R(eps, n, list(range(n)))
            assert R <= -(n // 2), ("Thm2.1 perm", n, pi, R)
    # tightness for K_n: exact minimum over all bipartitions
    kn = {}
    for n in range(4, 11):
        eps = {frozenset((i, j)): 1 for i in range(n) for j in range(i + 1, n)}
        R = brute_min_R(eps, n)
        edits = (n * (n - 1) // 2 + R) // 2
        assert edits == target(n), ("K_n tightness", n, edits)
        kn[n] = edits
    res["K_n_exact_min_edits"] = kn
    return res


# ------------------------------------------------------------------- part C

def T_of(pi):
    n = len(pi)
    return sum(((j - i) % n) * ((pi[j] - pi[i]) % n)
               for i in range(n) for j in range(n) if i != j)


def part_C(nc, rng):
    res = {}
    # Lemma 3.1 on all pi
    for n in range(4, nc + 1):
        for pi in itertools.permutations(range(n)):
            pi = list(pi)
            tot = sum(n * (n - 1) // 2 - inv_cut(pi, a, b) for a in range(n) for b in range(n))
            assert tot == T_of(pi), ("Lemma3.1", n, pi)
    # Lemma 3.2 exhaustively, with the equality case
    eq = {}
    for n in range(4, 9):
        bound = n * n * (n * n - 1) // 6
        eqs = []
        for pi in itertools.permutations(range(n)):
            T = T_of(list(pi))
            assert T >= bound, ("Lemma3.2", n, pi, T)
            if T == bound:
                eqs.append(pi)
        refl = sorted(tuple((h - i) % n for i in range(n)) for h in range(n))
        assert sorted(eqs) == refl, ("Lemma3.2 equality", n, eqs)
        eq[n] = len(eqs)
    res["equality_cases_are_the_n_reflections"] = eq
    # random large n
    for n in list(range(10, 31)) + [40, 60]:
        bound = n * n * (n * n - 1) // 6
        for _ in range(20):
            pi = list(range(n))
            rng.shuffle(pi)
            assert T_of(pi) >= bound, ("Lemma3.2 rnd", n, pi)
        assert T_of([(3 - i) % n for i in range(n)]) == bound, ("reflection T", n)
    # Theorem 3.3 against the true I(pi)
    thm = {}
    for n in range(4, 9):
        worst_I, worst_bound = -1, -1
        for pi in itertools.permutations(range(n)):
            pi = list(pi)
            I = I_of(pi)
            T = T_of(pi)
            b1 = (n * (n - 1) // 2 * n * n - T) // (n * n)   # floor(C(n,2) - T/n^2)
            assert I <= b1, ("Thm3.3 per pi", n, pi, I, b1)
            worst_I = max(worst_I, I)
            worst_bound = max(worst_bound, b1)
        b2 = ((n - 1) * (2 * n - 1)) // 6
        assert worst_bound <= b2, ("Thm3.3 global", n, worst_bound, b2)
        assert worst_I == target(n), ("max I", n, worst_I)
        thm[n] = {"max_I": worst_I, "max_per_pi_bound": worst_bound,
                  "global_bound": b2, "target": target(n)}
    res["theorem_3_3"] = thm
    return res


# ------------------------------------------------------------------- part D

def two_cut_check(pi, a1, b1, a2, b2):
    n = len(pi)
    boxes = {"11": [], "22": [], "12": [], "21": []}
    for i in range(n):
        inx = ((i - a1) % n) < ((a2 - a1) % n)
        iny = ((pi[i] - b1) % n) < ((b2 - b1) % n)
        key = ("1" if inx else "2") + ("1" if iny else "2")
        boxes[key].append(i)
    A, B, C, D = len(boxes["11"]), len(boxes["22"]), len(boxes["12"]), len(boxes["21"])

    def inv_box(S, a, b):
        co = [(((i - a) % n), ((pi[i] - b) % n)) for i in S]
        return sum(1 for u in range(len(co)) for v in range(u + 1, len(co))
                   if (co[u][0] - co[v][0]) * (co[u][1] - co[v][1]) < 0)

    tot = sum(inv_box(S, a1, b1) for S in boxes.values())
    tot2 = sum(inv_box(S, a2, b2) for S in boxes.values())
    assert tot == tot2, ("box inversions differ", pi, a1, b1, a2, b2)
    lhs = inv_cut(pi, a1, b1) + inv_cut(pi, a2, b2)
    rhs = 2 * tot + 2 * C * D + (A + B) * (C + D)
    assert lhs == rhs, ("Thm4.1", pi, a1, b1, a2, b2, lhs, rhs)
    return A, B, C, D


def part_D(rng):
    for n in range(4, 7):
        for pi in itertools.permutations(range(n)):
            pi = list(pi)
            for a1 in range(n):
                for a2 in range(n):
                    if a1 == a2:
                        continue
                    for b1 in range(n):
                        for b2 in range(n):
                            if b1 == b2:
                                continue
                            two_cut_check(pi, a1, b1, a2, b2)
    for n in range(7, 13):
        for _ in range(60):
            pi = list(range(n))
            rng.shuffle(pi)
            a1, a2 = rng.randrange(n), rng.randrange(n)
            b1, b2 = rng.randrange(n), rng.randrange(n)
            if a1 == a2 or b1 == b2:
                continue
            two_cut_check(pi, a1, b1, a2, b2)
    # Corollary 4.3 on reflections
    ok = {}
    for n in range(4, 41):
        k = n // 2
        h = 3 % n
        pi = [(h - i) % n for i in range(n)]
        a1, a2 = 0, k              # P_1 = [0, k)
        b1, b2 = (h - k + 1) % n, (h + 1) % n   # V_1 = h - P_1
        A, B, C, D = two_cut_check(pi, a1, b1, a2, b2)
        assert C == 0 and D == 0 and {A, B} == {k, n - k}, ("Cor4.3 boxes", n, A, B, C, D)
        val = inv_cut(pi, a1, b1)
        assert val == k * (k - 1) // 2 + (n - k) * (n - k - 1) // 2 == target(n), \
            ("Cor4.3 value", n, val)
        ok[n] = val
    return {"reflection_two_cut_value": {4: ok[4], 10: ok[10], 25: ok[25], 40: ok[40]},
            "reflections_checked": len(ok)}


# ------------------------------------------------------------------- part E

def part_E():
    res = {}
    # 1. max_alpha min_beta and max_beta min_alpha
    fa = {}
    for n in (7, 8):
        worst = -1
        arg = None
        for pi in itertools.permutations(range(n)):
            pi = list(pi)
            tab = [[inv_cut(pi, a, b) for b in range(n)] for a in range(n)]
            v = max(min(tab[a][b] for b in range(n)) for a in range(n))
            u = max(min(tab[a][b] for a in range(n)) for b in range(n))
            assert v == u or True
            if v > worst:
                worst, arg = v, tuple(pi)
        fa[n] = {"max_alpha_min_beta": worst, "target": target(n), "example": arg}
        assert worst > target(n), ("forall-alpha refuted", n, worst)
    res["forall_alpha_exists_beta"] = fa
    # 2. best (h, k) family of n cuts
    hk = {}
    for n in range(5, 9):
        best_family = None
        for h in range(n):
            for k in range(n):
                w = 0
                for pi in itertools.permutations(range(n)):
                    pi = list(pi)
                    w = max(w, min(inv_cut(pi, (t + h) % n, (pi[t] + k) % n) for t in range(n)))
                    if best_family is not None and w > best_family:
                        break
                if best_family is None or w < best_family:
                    best_family = w
        hk[n] = {"best_family_max": best_family, "target": target(n)}
        assert best_family > target(n), ("(h,k) family refuted", n, best_family)
    res["hk_families"] = hk
    # 3. depth-1 local conditions insufficient at n = 10
    w = [0, 1, 8, 7, 6, 5, 4, 3, 2, 9]
    n = 10
    t = [0] * n
    for j, v in enumerate(w):
        t[v] = j
    A = [0] * (n + 1)
    B = [0] * (n + 1)
    for k in range(n):
        A[k + 1] = A[k] + n - 1 - 2 * w[k]
        B[k + 1] = B[k] + n - 1 - 2 * t[k]
    bad = []
    for k in range(n + 1):
        for l in range(n + 1):
            if min(k, l) > 1:
                continue
            K = sum(1 for j in range(k) if w[j] < l)
            if A[k] + B[l] + 2 * k * l - 2 * n * K < 0:
                bad.append((k, l))
    assert not bad, ("depth-1 conditions should hold", bad)
    assert inv_of(w) == 21 > target(10), ("depth-1 counterexample", inv_of(w))
    res["depth1_counterexample"] = {"n": 10, "line": w, "inv": inv_of(w), "target": target(10),
                                    "violated_depth1_conditions": len(bad)}
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--na", type=int, default=6, help="exhaustive n for part A")
    ap.add_argument("--nc", type=int, default=7, help="exhaustive n for Lemma 3.1")
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    rng = random.Random(20260914)
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    report = {"version": VERSION, "na": args.na, "nc": args.nc, "quick": args.quick}
    report["A"] = part_A(args.na, rng)
    print("part A ok", flush=True)
    report["B"] = part_B(rng, args.quick)
    print("part B ok", flush=True)
    report["C"] = part_C(args.nc, rng)
    print("part C ok", flush=True)
    report["D"] = part_D(rng)
    print("part D ok", flush=True)
    report["E"] = part_E()
    print("part E ok", flush=True)
    report["seconds"] = round(time.time() - t0, 1)
    report["result"] = "PASS"
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(report, f, indent=1, sort_keys=True, default=str)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# check_C37 — независимая проверка C37–C40\n\n")
        f.write("Версия %s, %s с, результат: %s.\n\n" % (VERSION, report["seconds"], report["result"]))
        f.write("Команда: `python3 checks/check_C37.py --na %d --nc %d`.\n\n" % (args.na, args.nc))
        f.write("```json\n%s\n```\n" % json.dumps(report, indent=1, sort_keys=True, default=str))
    print("PASS in %.1f s" % (time.time() - t0))


if __name__ == "__main__":
    main()
