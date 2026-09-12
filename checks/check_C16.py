"""Independent finite checks of N1 (docs/incoming/LRX_STRICT_PROOF_RU-1.md).

Written from the text of N1 after the logical audit (docs/notes/
strict_proof_audit.md). This is NOT the author's verify_lrx.py (not available):
loads are counted by explicit edge enumeration, contraction is checked by the
actual presence of edge {0,1} in each arc, and every word is applied letter by
letter with the reference moves of oracle/moves.py (oracle-1.0).

Each check names the N1 equation it tests. Passing finite checks do not prove
the general theorem; they test the implementation of the construction and hunt
for errors in the lemmas on small inputs (AGENTS.md rules 5, 9).

Usage: python checks/check_C16.py [--cycle-max-n 9] [--perm-max-n 7]
       [--scalar-max-n 100000] [--quick]
Output: data/runs/check_C16/report.json and summary.md
"""

import argparse
import itertools
import json
import os
import random
import sys
import time
from fractions import Fraction
from math import isqrt

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
sys.path.insert(0, os.path.join(ROOT, "exact"))
sys.path.insert(0, os.path.join(ROOT, "bounds"))

from moves import apply_word, freely_reduce, identity, sigma, rev, delta, CORE_VERSION  # noqa: E402
from bfs import rank_perm, factorials  # noqa: E402
import known  # noqa: E402

OUT_DIR = os.path.join(ROOT, "data", "runs", "check_C16")


class CheckFailure(Exception):
    pass


def check(cond, tag, info=None):
    if not cond:
        raise CheckFailure(f"{tag}: {info!r}")


# ---------------------------------------------------------------- N1 §1

def signed_step(a, b, n):
    """d = b - a mod n with -n/2 < d <= n/2; antipodal step positive."""
    d = (b - a) % n
    if 2 * d > n:
        d -= n
    return d


def arc_edges(a, d, n):
    """Edge indices v (edge {v, v+1}) traversed by the arc from a with signed step d."""
    if d > 0:
        return [(a + i) % n for i in range(d)]
    return [(a - 1 - i) % n for i in range(-d)]


def cycle_data(cycle, n):
    """Loads, F, M, t, E, S for an oriented nontrivial cycle (list of positions)."""
    k = len(cycle)
    steps = [signed_step(cycle[j], cycle[(j + 1) % k], n) for j in range(k)]
    loads = [0] * n
    for j in range(k):
        for v in arc_edges(cycle[j], steps[j], n):
            loads[v] += 1
    F = sum(abs(d) for d in steps)
    M = max(loads)
    t = sum(1 for j in range(k) if steps[j - 1] < 0 < steps[j])
    t_pm = sum(1 for j in range(k) if steps[j - 1] > 0 > steps[j])
    E = k - 2 * t
    S = 2 * M + E - 1
    return {"k": k, "steps": steps, "loads": loads, "F": F, "M": M, "t": t,
            "t_pm": t_pm, "E": E, "S": S}


def support_connected(loads, n):
    pos = [loads[v] > 0 for v in range(n)]
    if all(pos):
        return True
    changes = sum(1 for v in range(n) if pos[v] != pos[v - 1])
    return any(pos) and changes == 2


def check_section1(cycle, n, cd):
    """(1.1) jump rule, sum of loads = F, F <= nM, t(+-) = t(-+), E >= 0, support connected."""
    k = cd["k"]
    steps, loads = cd["steps"], cd["loads"]
    check(all(d != 0 for d in steps), "steps nonzero", cycle)
    check(sum(loads) == cd["F"], "sum loads = F", cycle)
    check(cd["F"] <= n * cd["M"], "F <= nM", cycle)
    check(cd["t_pm"] == cd["t"], "t(+->-) = t(-->+)", cycle)
    check(cd["E"] >= 0, "E >= 0", cycle)
    pos = {cycle[j]: j for j in range(k)}
    for v in range(n):
        jump = loads[v] - loads[v - 1]
        if v not in pos:
            check(jump == 0, "(1.1) jump outside C", (cycle, v))
        else:
            j = pos[v]
            din, dout = steps[j - 1], steps[j]
            if din < 0 < dout:
                exp = 2
            elif din > 0 > dout:
                exp = -2
            else:
                exp = 0
            check(jump == exp, "(1.1) jump at cycle vertex", (cycle, v, jump, exp))
    check(support_connected(loads, n), "positive support connected", cycle)


# ---------------------------------------------------------------- N1 §2

def choose_carrier(cycle, n, cd):
    """§2.1: index j0 with a_{j0} = b such that edge {b, b+1} has load M; mixed case
    b is the start of a maximal run of load M and a -> + vertex."""
    steps, loads, M = cd["steps"], cd["loads"], cd["M"]
    signs = {d > 0 for d in steps}
    if len(signs) == 1:
        check(all(l == M for l in loads), "unidirectional: constant loads", cycle)
        return 0, "unidirectional"
    starts = [v for v in range(n) if loads[v] == M and loads[v - 1] < M]
    check(len(starts) >= 1, "mixed: a maximal M-run exists", cycle)
    b = starts[0]
    check(b in cycle, "carrier in cycle", (cycle, b))
    j0 = cycle.index(b)
    check(steps[j0 - 1] < 0 < steps[j0], "(2.1) carrier is a -->+ vertex", (cycle, b))
    return j0, "mixed"


def positional_perm_of_word(word, n):
    """sigma with apply_word(p, word) == p o sigma, i.e. sigma[i] = position read
    into slot i; obtained by applying the word to the identity array."""
    return apply_word(identity(n), word)


def local_word(cycle, n):
    """§2.2–2.4: build W_C in carrier coordinates (carrier at 0), return details.

    Returns dict with e_j, raw word, cancelled word, expected positional
    permutation (0 v_1 ... v_{k-1})^{-1} and checks (2.2)–(2.6).
    """
    cd = cycle_data(cycle, n)
    j0, case = choose_carrier(cycle, n, cd)
    k = cd["k"]
    cyc = cycle[j0:] + cycle[:j0]
    b = cyc[0]
    v = [(a - b) % n for a in cyc]           # carrier coordinates, v_0 = 0
    steps = [signed_step(v[j], v[(j + 1) % k], n) for j in range(k)]
    check(steps == cd["steps"][j0:] + cd["steps"][:j0], "steps invariant under shift", cycle)
    e = []
    for j in range(k):
        edges = arc_edges(v[j], steps[j], n)
        passes = 0 in edges                    # edge {0,1}
        e.append(steps[j] - (1 if steps[j] > 0 else -1) if passes else steps[j])
    F, M, t, E, S = cd["F"], cd["M"], cd["t"], cd["E"], cd["S"]
    check(sum(abs(x) for x in e) == F - M, "(2.2) sum|e| = F - M", (cycle, e))
    check(sum(e) % (n - 1) == 0, "(2.3) sum e = 0 mod n-1", (cycle, e))
    for j in range(1, k - 1):
        check(e[j] != 0 and (e[j] > 0) == (steps[j] > 0), "inner e_j nonzero, sign of d_j", (cycle, j, e))
    for j in range(k):
        if e[j] == 0:
            check((j == 0 and steps[0] > 0) or (j == k - 1 and steps[k - 1] < 0),
                  "only first positive / last negative step may vanish", (cycle, j))
    # projection identity u^{z_j}(1) = v_j on the ring 1..n-1
    z = 0
    for j in range(1, k):
        z += e[j - 1]
        check((z % (n - 1)) + 1 == v[j] if v[j] != 0 else False, "u^{z_j}(1) = v_j", (cycle, j, z, v))
    # raw word (2.4): blocks U^{e_j}, X between blocks
    blocks = []
    for x in e:
        blocks.append(list("XL" * x) if x > 0 else list("RX" * (-x)))
    tokens = []       # (letter, block_index, position_in_block) ; inserted X marked block -1
    for j in range(k):
        tokens.extend((ch, j, i) for i, ch in enumerate(blocks[j]))
        if j < k - 1:
            tokens.append(("X", -1, j + 1))    # inserted X for vertex a_{j+1}
    raw = "".join(tok[0] for tok in tokens)
    check(len(raw) == 2 * (F - M) + k - 1, "raw length 2(F-M)+k-1", cycle)
    check(raw.count("X") == (F - M) + k - 1, "raw N_X", cycle)
    # cancellations per §2.4
    kill = set()
    ncancel = 0
    for j in range(1, k):                      # inserted X for vertex a_j
        ins_idx = next(i for i, tok in enumerate(tokens) if tok[1] == -1 and tok[2] == j)
        dout, din = steps[j], steps[j - 1]
        if dout > 0:
            check(e[j] > 0, "next block nonempty", (cycle, j))
            nxt = ins_idx + 1
            check(tokens[nxt][0] == "X" and tokens[nxt][1] == j and tokens[nxt][2] == 0, "first X of next block", cycle)
            pair = (ins_idx, nxt)
        elif din < 0:
            check(e[j - 1] < 0, "previous block nonempty", (cycle, j))
            prv = ins_idx - 1
            check(tokens[prv][0] == "X" and tokens[prv][1] == j - 1
                  and tokens[prv][2] == len(blocks[j - 1]) - 1, "last X of previous block", cycle)
            pair = (prv, ins_idx)
        else:
            check(din > 0 > dout, "+->- vertex has no cancellation", cycle)
            continue
        check(not (set(pair) & kill), "cancellations disjoint", (cycle, pair))
        kill.update(pair)
        ncancel += 1
    check(ncancel == k - 1 - t, "number of cancellations = k-1-t", (cycle, ncancel, k, t))
    word = "".join(tok[0] for i, tok in enumerate(tokens) if i not in kill)
    check(len(word) == 2 * F - S, "(2.6) |W_C| = 2F - S after chosen cancellations", (cycle, len(word), 2 * F - S))
    # action (2.5): positional permutation equals (0 v_1 ... v_{k-1})^{-1}
    cyc_perm = list(range(n))
    for j in range(k):
        cyc_perm[v[j]] = v[(j + 1) % k]        # C'(v_j) = v_{j+1}
    inv = [0] * n
    for i in range(n):
        inv[cyc_perm[i]] = i
    expected = tuple(inv)                      # id o C'^{-1}
    check(positional_perm_of_word(raw, n) == expected, "(2.5) raw word action", (cycle, raw))
    check(positional_perm_of_word(word, n) == expected, "(2.5) cancelled word action", (cycle, word))
    reduced = freely_reduce(word)
    check(len(reduced) <= len(word), "free reduction", cycle)
    return {"case": case, "carrier": b, "v": v, "e": e, "raw": raw, "word": word,
            "reduced_len": len(reduced), "cd": cd}


# ---------------------------------------------------------------- N1 §7

def check_cycle_bound(cycle, n, cd):
    """(7.1)–(7.2) and the small-load classification of §7.1."""
    k, F, M, t, E, S = cd["k"], cd["F"], cd["M"], cd["t"], cd["E"], cd["S"]
    steps = cd["steps"]
    antipodal_transp = (k == 2 and n % 2 == 0 and steps[0] == n // 2 and steps[1] == n // 2)
    if M == 1:
        check(t == 0 and E == k, "M=1: no turns", cycle)
        check(k >= 3 or antipodal_transp, "M=1,k=2 only antipodal transposition", cycle)
    if M <= 2 and t > 0:
        check(M == 2 and set(cd["loads"]) == {0, 2}, "M<=2 with turn: loads in {0,2}", cycle)
        check(t == 1 and k == E + 2, "(7.2) t=1, k=E+2", cycle)
        if E == 0:
            check(k == 2 and F <= n and not antipodal_transp, "(7.2) M=2,E=0: non-antipodal transposition", cycle)
        if E == 1:
            check(k == 3 and F <= n and sum(steps) == 0, "(7.2) M=2,E=1: triangle, F<=n", cycle)
    rhs = Fraction(4, 3) * M + Fraction(F, 3 * n) + Fraction(7, 9) * E - (Fraction(2, 9) if antipodal_transp else 0)
    check(S >= rhs, "(7.1) S >= 4M/3 + F/(3n) + 7E/9 - 2/9*1", (cycle, S, rhs))
    if antipodal_transp:
        check(S == rhs, "(7.1) equality on antipodal transposition", cycle)
    return antipodal_transp


def all_oriented_cycles(n):
    """All nontrivial oriented cycles on Z_n, one representative per cyclic rotation
    (the smallest element first)."""
    for k in range(2, n + 1):
        for subset in itertools.combinations(range(n), k):
            first = subset[0]
            for rest in itertools.permutations(subset[1:]):
                yield [first] + list(rest)


def run_cycles(max_n):
    counts = {}
    for n in range(4, max_n + 1):
        cnt = 0
        antipodal = 0
        for cyc in all_oriented_cycles(n):
            cd = cycle_data(cyc, n)
            check_section1(cyc, n, cd)
            local_word(cyc, n)
            if check_cycle_bound(cyc, n, cd):
                antipodal += 1
            cnt += 1
        counts[n] = {"oriented_cycles": cnt, "antipodal_transpositions": antipodal}
    return counts


# ---------------------------------------------------------------- N1 §§3–6, 8

def cycles_of(f):
    n = len(f)
    seen = [False] * n
    out = []
    for i in range(n):
        if seen[i] or f[i] == i:
            seen[i] = True
            continue
        cyc = []
        j = i
        while not seen[j]:
            seen[j] = True
            cyc.append(j)
            j = f[j]
        out.append(cyc)
    return out


def H_of(c, n):
    return n if c % n == 0 else n + delta(0, c, n) - 2


def external_route(c, n):
    """§3: list of external positions visited (starting at 0) and the letters."""
    c %= n
    if c == 0:
        return [0] + [(i + 1) % n for i in range(n)], "L" * n
    d = delta(0, c, n)
    dirn = 1 if c <= n // 2 else -1
    pos, letters = [0], []
    cur = 0
    for _ in range(d - 1):
        cur = (cur + dirn) % n; pos.append(cur); letters.append("L" if dirn > 0 else "R")
    for _ in range(d - 1):
        cur = (cur - dirn) % n; pos.append(cur); letters.append("L" if dirn < 0 else "R")
    for _ in range(n - d):
        cur = (cur - dirn) % n; pos.append(cur); letters.append("L" if dirn < 0 else "R")
    check(cur == c, "route ends at c", (c, n))
    check(len(letters) == H_of(c, n), "route length = H(c)", (c, n))
    check(set(pos) == set(range(n)), "route visits all vertices", (c, n))
    return pos, "".join(letters)


def shift_stats(pi, c, n, build_word=True):
    """All per-shift statistics of rule 15 (AGENTS.md) and the full sorting word."""
    f = [(pi[i] + c) % n for i in range(n)]
    cycs = cycles_of(f)
    F = S = sumM = E = theta = 0
    all_loads = [0] * n
    words = {}
    local_total = 0
    for cyc in cycs:
        cd = cycle_data(cyc, n)
        lw = local_word(cyc, n)
        F += cd["F"]; S += cd["S"]; sumM += cd["M"]; E += cd["E"]
        for v in range(n):
            all_loads[v] += cd["loads"][v]
        k, steps = cd["k"], cd["steps"]
        if k == 2 and n % 2 == 0 and steps[0] == n // 2 and steps[1] == n // 2:
            theta += 1
        check(lw["carrier"] not in words, "distinct carriers", (pi, c))
        words[lw["carrier"]] = lw["word"]
        local_total += len(lw["word"])
    check(local_total == 2 * F - S, "sum of local words = 2F_c - S_c", (pi, c))
    M_all = max(all_loads)
    check(M_all <= sumM, "M_c^all <= sum_C M(C)", (pi, c))
    check(F <= n * sumM, "F_c <= n sum M", (pi, c))
    res = {"c": c, "F_c": F, "S_c": S, "sumM": sumM, "E_c": E, "theta_c": theta,
           "H": H_of(c, n), "M_all": M_all, "all_loads": all_loads,
           "n_cycles": len(cycs), "bound": 2 * F - S + H_of(c, n)}
    if build_word:
        pos, letters = external_route(c, n)
        full = []
        done = set()
        for idx, b in enumerate(pos):
            if b in words and b not in done:
                full.append(words[b]); done.add(b)
            if idx < len(letters):
                full.append(letters[idx])
        check(done == set(words), "all carriers visited", (pi, c))
        word = "".join(full)
        check(len(word) == 2 * F - S + H_of(c, n), "(3.3) full word length", (pi, c, len(word)))
        check(apply_word(tuple(pi), word) == identity(n), "(3.2) full word sorts pi", (pi, c, word))
        red = freely_reduce(word)
        nx = word.count("X")
        res.update({"len": len(word), "N_X": nx, "N_rot": len(word) - nx,
                    "reduced_len": len(red), "reduced_N_X": red.count("X"),
                    "route_len": len(letters), "word": word})
    return res


def half_window(a, n):
    return {(a + i) % n for i in range(1, n // 2 + 1)}


def perm_stats(pi, build_words=True, shifts=None, exact_d=None):
    """Checks of §§3–6, 8 for one permutation; returns statistics."""
    n = len(pi)
    P = known.P(n)
    tau = [(-pi[i]) % n for i in range(n)]
    shifts_all = list(range(n))
    per_c = {}
    for c in shifts_all:
        per_c[c] = shift_stats(pi, c, n, build_word=build_words and (shifts is None or c in shifts))
    sumF = sum(per_c[c]["F_c"] for c in shifts_all)
    sumH = sum(per_c[c]["H"] for c in shifts_all)
    check(sumF == n * P, "bar F = P", pi)
    check(Fraction(sumH, n) == n - 2 + Fraction(P + 2, n), "bar H = n-2+(P+2)/n", pi)
    s = Fraction(sum(per_c[c]["S_c"] for c in shifts_all), n)
    mu = Fraction(sum(per_c[c]["sumM"] for c in shifts_all), n)
    Ebar = Fraction(sum(per_c[c]["E_c"] for c in shifts_all), n)
    thetabar = Fraction(sum(per_c[c]["theta_c"] for c in shifts_all), n)
    calE = sum(per_c[c]["E_c"] for c in shifts_all)
    expr34 = 2 * P + n - 2 + Fraction(P + 2, n) - s
    min_bound = min(per_c[c]["bound"] for c in shifts_all)
    check(min_bound <= expr34, "(3.4) min_c bound <= averaged expression", pi)
    # §4: e_i, V, B_i, identity (4.3)-(4.5)
    windows = [half_window(i, n) for i in range(n)]
    e = []
    Bsets = []
    for i in range(n):
        tI = {tau[j] for j in windows[i]}
        e.append(len(tI ^ windows[tau[i]]))
        Bi = {j for j in range(n) if (j in windows[i]) != (tau[j] in windows[tau[i]])}
        check(len(Bi) == e[i], "(4.1) |B_i| = e_i", (pi, i))
        Bsets.append(Bi)
    V = sum(e)
    x = Fraction(V, n)
    # (4.3): calE = #{(i,j): i != j, sgn(j-i) = sgn(pi(j)-pi(i))}
    cnt43 = sum(1 for i in range(n) for j in range(n) if i != j
                and (signed_step(i, j, n) > 0) == (signed_step(pi[i], pi[j], n) > 0))
    check(cnt43 == calE, "(4.3) calE by sign count", pi)
    B = 0
    if n % 2 == 0:
        B = sum(1 for i in range(n) for j in range(i + 1, n)
                if (j - i) % n == n // 2 and (pi[j] - pi[i]) % n == n // 2)
    check(calE == V + 2 * B, "(4.4) calE = V + 2B", (pi, calE, V, B))
    check(B == sum(per_c[c]["theta_c"] for c in shifts_all), "(4.5) B = sum theta_c", pi)
    check(x == Ebar - 2 * thetabar, "(4.5) x = bar E - 2 bar theta", pi)
    # §5: metric identity, (5.2), (5.3), exceptional pairs, anchor
    u = [(tau[i] - i) % n for i in range(n)]
    for a in range(n):
        for b in range(n):
            check(len(windows[a] ^ windows[b]) == 2 * delta(a, b, n), "|I_a ^ I_b| = 2 delta", (n, a, b))
            if a == b:
                continue
            dab, dtt = delta(a, b, n), delta(tau[a], tau[b], n)
            check(2 * abs(dab - dtt) <= e[a] + e[b], "(5.2)", (pi, a, b))
            s1, s2 = signed_step(a, b, n), signed_step(tau[a], tau[b], n)
            exceptional = ((s1 > 0) != (s2 > 0)) and 2 * abs(s1) != n and 2 * abs(s2) != n
            if not exceptional:
                check(delta(u[a], u[b], n) == abs(dab - dtt), "(5.3) equality for non-exceptional pair", (pi, a, b))
            else:
                check(b in Bsets[a] and a in Bsets[b], "exceptional pair mutual in B", (pi, a, b))
    def T(h):
        return sum(delta((pi[i] + i) % n, h, n) for i in range(n))
    m = min(e)
    Q = V - n * m
    anchor_ok = []
    for a in [i for i in range(n) if e[i] == m]:
        h = (a + pi[a]) % n
        Th = T(h)
        check(Th == sum(delta(u[i], u[a], n) for i in range(n)), "(5.4)", (pi, a))
        if m == 0:
            check(2 * Th <= V, "(5.5) T_h <= V/2", (pi, a, Th, V))
        else:
            for j in range(n):
                if j == a:
                    continue
                s1, s2 = signed_step(a, j, n), signed_step(tau[a], tau[j], n)
                exceptional = ((s1 > 0) != (s2 > 0)) and 2 * abs(s1) != n and 2 * abs(s2) != n
                if exceptional:
                    rho_j = e[j] - m
                    check(delta(u[j], u[a], n) <= 2 * m + Fraction(rho_j, 2) + Fraction(Q, 2 * m), "(5.6)", (pi, a, j))
            check(Th <= V + m * (m - 1), "(5.8) T_h <= V + m(m-1)", (pi, a, Th))
        check(Th <= (n - 1) * x + x * x, "(5.1) T_h <= (1-1/n)V + V^2/n^2", (pi, a, Th))
        anchor_ok.append((a, h, Th))
    a0, h0, Th0 = anchor_ok[0]
    # §6: (6.1) for all h, with (6.2)-(6.3) through the explicit window a=(h+c)/2
    twoSumM = 2 * sum(per_c[c]["M_all"] for c in shifts_all)
    inv2 = pow(2, -1, n) if n % 2 else None
    for h in range(n):
        Th = T(h)
        check(twoSumM >= n * (n - 1) - 2 * Th, "(6.1)", (pi, h))
        total_D = 0
        for c in shifts_all:
            a = ((h + c) * inv2) % n if n % 2 else ((h + c) // 2) % n
            W = windows[a]
            D = sum(1 for i in range(n) if (i in W) != (((pi[i] + c) % n) in W))
            loads = per_c[c]["all_loads"]
            check(D <= loads[a] + loads[(a + n // 2) % n] <= 2 * per_c[c]["M_all"], "(6.2)", (pi, h, c))
            total_D += D
        check(total_D >= n * (n - 1) - 2 * Th, "(6.3) summed separations", (pi, h))
    check(mu >= Fraction(P, n), "(6.4a) mu >= P/n", pi)
    check(mu >= Fraction(n - 1, 2) - Fraction(Th0, n), "(6.4b) mu >= (n-1)/2 - T_h/n", pi)
    # §7 (7.3) and §8
    check(s >= Fraction(4, 3) * mu + Fraction(P, 3 * n) + Fraction(7, 9) * x + Fraction(4, 3) * thetabar, "(7.3)", pi)
    A = Fraction(5 * P, 3 * n) + Fraction(7, 9) * x
    Bx = Fraction(2 * (n - 1), 3) + Fraction(P, 3 * n) - (Fraction(5, 9) - Fraction(4, 3 * n)) * x - Fraction(4, 3 * n) * x * x
    check(s >= A, "(8.1)", pi)
    check(s >= Bx, "(8.2)", pi)
    # (8.3): s >= 5P/(3n) + 7 r_n / 9  <=>  y := (9/7)(s - 5P/(3n)) >= r_n
    y = Fraction(9, 7) * (s - Fraction(5 * P, 3 * n))
    z = 3 * n * n - 4 * n + 1 - 4 * P
    lhs = 2 * y + (n - 1)
    check(lhs >= 0 and lhs * lhs >= z, "(8.3) s >= 5P/(3n) + 7 r_n/9", (pi, s))
    U = known.strict_upper_bound(n)
    # final: min_c actual length <= floor(expr34) <= U_n
    out = {"pi": list(pi), "n": n, "V": V, "x": str(x), "m": m, "Q": Q, "calE": calE, "B": B,
           "s": str(s), "mu": str(mu), "barE": str(Ebar), "bar_theta": str(thetabar),
           "expr_3_4": str(expr34), "min_bound_c": min_bound, "U_n": U,
           "anchor": {"a": a0, "h": h0, "T_h": Th0}, "T_min": min(T(h) for h in range(n))}
    check(min_bound <= expr34 < U + 1, "min_c bound <= (3.4) < U_n + 1", (pi, min_bound, expr34, U))
    built = [c for c in shifts_all if "len" in per_c[c]]
    if built:
        lens = {c: per_c[c]["len"] for c in built}
        red = {c: per_c[c]["reduced_len"] for c in built}
        cmin = min(lens, key=lens.get)
        out.update({"shifts_built": len(built), "min_len": lens[cmin], "argmin_c": cmin,
                    "min_reduced_len": min(red.values()),
                    "N_X_at_argmin": per_c[cmin]["N_X"], "N_rot_at_argmin": per_c[cmin]["N_rot"],
                    "mean_len": str(Fraction(sum(lens.values()), len(lens)))})
        for c in built:
            check(lens[c] == per_c[c]["bound"], "len == 2F_c - S_c + H(c)", (pi, c))
    if exact_d is not None:
        out["d_exact"] = exact_d
        check(out.get("min_len", min_bound) >= exact_d, "word length >= exact distance", pi)
        out["loss_min_len"] = out["min_len"] - exact_d
        out["loss_min_reduced"] = out["min_reduced_len"] - exact_d
    out["per_c"] = {c: {k: v for k, v in per_c[c].items() if k not in ("all_loads", "word")} for c in shifts_all}
    return out


def load_table(n):
    path = os.path.join(ROOT, "data", "tables", f"dist_n{n}.bin")
    with open(path, "rb") as f:
        return f.read()


def run_perms_exhaustive(max_n):
    summary = {}
    for n in range(4, max_n + 1):
        table = load_table(n)
        fact = factorials(n)
        losses = {}
        losses_red = {}
        worst = []
        cnt = 0
        max_expr_minus_B = None
        for pi in itertools.permutations(range(n)):
            d = table[rank_perm(pi, fact)]
            st = perm_stats(pi, build_words=True, exact_d=d)
            cnt += 1
            losses[st["loss_min_len"]] = losses.get(st["loss_min_len"], 0) + 1
            losses_red[st["loss_min_reduced"]] = losses_red.get(st["loss_min_reduced"], 0) + 1
            worst.append((st["loss_min_len"], list(pi), st["min_len"], d))
            gap = Fraction(st["expr_3_4"]) - known.target_diameter(n)
            if max_expr_minus_B is None or gap > max_expr_minus_B:
                max_expr_minus_B = gap
        worst.sort(reverse=True)
        summary[n] = {"permutations": cnt,
                      "loss_min_len_hist": dict(sorted(losses.items())),
                      "loss_min_reduced_hist": dict(sorted(losses_red.items())),
                      "max_loss_min_len": worst[0][0],
                      "worst_examples": worst[:5],
                      "max_expr_3_4_minus_B_n": str(max_expr_minus_B),
                      "U_n": known.strict_upper_bound(n), "D_n": max(table)}
    return summary


def reflections(n):
    return [tuple((h - i) % n for i in range(n)) for h in range(n)]


def run_reflections(max_n):
    """C24: on reflections s = 3(n-1)/2 exactly and (3.4) = B_n + P/n - 1/2 + 2/n (even),
    B_n + P/n - 1 + 2/n (odd). Also rule-15 statistics for sigma_n and all its rotations."""
    out = {}
    for n in range(4, max_n + 1):
        Bn, P = known.target_diameter(n), known.P(n)
        expected = Bn + Fraction(P, n) + Fraction(2, n) - (Fraction(1, 2) if n % 2 == 0 else 1)
        rows = []
        for pi in reflections(n):
            st = perm_stats(pi, build_words=True)
            check(Fraction(st["s"]) == Fraction(3 * (n - 1), 2), "C24: s = 3(n-1)/2 on reflections", pi)
            check(Fraction(st["expr_3_4"]) == expected, "C24: (3.4) value on reflections", pi)
            rows.append({"pi": list(pi), "min_len": st["min_len"], "argmin_c": st["argmin_c"],
                         "N_X": st["N_X_at_argmin"], "N_rot": st["N_rot_at_argmin"],
                         "min_reduced_len": st["min_reduced_len"], "mean_len": st["mean_len"],
                         "T_min": st["T_min"], "V": st["V"]})
        sig = perm_stats(sigma(n), build_words=True)
        out[n] = {"B_n": Bn, "expr_3_4_on_reflections": str(expected), "rows": rows,
                  "sigma_n_per_c": sig["per_c"], "sigma_n_min_len": sig["min_len"],
                  "sigma_n_min_reduced_len": sig["min_reduced_len"]}
    return out


def structured_samples(n, rng, k_random):
    """Rule 9 sample families."""
    out = {"sigma_n": sigma(n), "rev_n": rev(n)}
    out["sigma_rot_L1"] = apply_word(sigma(n), "L")
    out["sigma_rot_L2"] = apply_word(sigma(n), "LL")
    cyc = tuple((i + 1) % n for i in range(n))
    out["long_cycle"] = cyc
    out["long_cycle_inv"] = tuple((i - 1) % n for i in range(n))
    out["short_cycles_3"] = tuple((i // 3) * 3 + (i % 3 + 1) % 3 if i // 3 < n // 3 else i for i in range(n))
    tp = list(range(n))
    for i in range(0, n - 1, 2):
        tp[i], tp[i + 1] = tp[i + 1], tp[i]
    out["adjacent_transpositions"] = tuple(tp)
    if n % 2 == 0:
        out["antipodal_transpositions"] = tuple((i + n // 2) % n for i in range(n))
    p = list(sigma(n))
    i, j = rng.randrange(n), rng.randrange(n)
    p[i], p[j] = p[j], p[i]
    out["sigma_perturbed"] = tuple(p)
    for r in range(k_random):
        q = list(range(n)); rng.shuffle(q)
        out[f"random_{r}"] = tuple(q)
    return out


def run_samples(spec, seed):
    rng = random.Random(seed)
    out = {}
    for n, k_random, shifts in spec:
        rows = {}
        for name, pi in structured_samples(n, rng, k_random).items():
            sh = None if shifts is None else set(range(n))  # None -> all shifts built
            st = perm_stats(pi, build_words=True, shifts=sh)
            rows[name] = {k: v for k, v in st.items() if k != "per_c"}
            rows[name]["pi"] = list(pi)
        out[n] = {"cases": len(rows), "seed": seed, "rows": rows}
    return out


# ---------------------------------------------------------------- N1 §10.1 (C23)

def run_family_C23(bs):
    """N1 §10.1 family for n = 8b: tau(bk+j) = b(3k mod 8) + j, pi = -tau.
    Claimed: e_i = 2b for all i, V = n^2/4, T_min = n^2/4, B = n/2, calE = V + n."""
    out = {}
    for b in bs:
        n = 8 * b
        tau = [0] * n
        for k in range(8):
            for j in range(b):
                tau[b * k + j] = b * ((3 * k) % 8) + j
        pi = tuple((-tau[i]) % n for i in range(n))
        check(sorted(pi) == list(range(n)), "C23 family is a permutation", n)
        st = perm_stats(pi, build_words=(n <= 16))
        windows = [half_window(i, n) for i in range(n)]
        e = [len({tau[j] for j in windows[i]} ^ windows[tau[i]]) for i in range(n)]
        check(all(v == 2 * b for v in e), "C23: e_i = 2b", (n, e))
        check(st["V"] == n * n // 4, "C23: V = n^2/4", (n, st["V"]))
        check(st["T_min"] == n * n // 4, "C23: T_min = n^2/4", (n, st["T_min"]))
        check(st["B"] == n // 2, "C23: B = n/2", (n, st["B"]))
        check(st["calE"] == st["V"] + n, "C23: calE = V + n", n)
        check(Fraction(st["T_min"], st["calE"]) == Fraction(n, n + 4), "C23: T_min/calE = n/(n+4)", n)
        out[n] = {"b": b, "V": st["V"], "T_min": st["T_min"], "B": st["B"], "calE": st["calE"],
                  "s": st["s"], "expr_3_4": st["expr_3_4"], "U_n": st["U_n"],
                  "min_len": st.get("min_len"), "d_note": "exact d unknown for n > 12"}
    return out


# ---------------------------------------------------------------- N1 §9

def run_scalar(max_n):
    """Exact integer checks of (T) -> (T') rounding for 4 <= n <= max_n, plus the
    constants of §9."""
    # q = (sqrt2-1)/2 < 5/24  <=> sqrt2 < 17/12 <=> 2 < 289/144
    check(2 * 144 < 289, "q < 5/24")
    # -5/2 + 35/216 + 13/42 = -3067/1512 < -2 ; -2 + 35/216 + 1/2 + 7/144 = -557/432 < -1
    check(Fraction(-5, 2) + Fraction(35, 216) + Fraction(13, 42) == Fraction(-3067, 1512) < -2, "odd constant")
    check(-2 + Fraction(35, 216) + Fraction(1, 2) + Fraction(7, 144) == Fraction(-557, 432) < -1, "even constant")
    # eps_4 = 1/(2(3 sqrt2 + sqrt17)) < 1/16  <=> 3 sqrt2 + sqrt17 > 8 ; 3 sqrt2 > 4 (18 > 16)
    check(18 > 16, "eps_4 < 1/16")
    # n = 5: 727/45 - 14 sqrt2/9 < 14  <=> (97/45)^2 < 392/81
    check(Fraction(97, 45) ** 2 < Fraction(392, 81), "R_5 < 14")
    check(known.strict_upper_bound(5) == 13 == known.strict_upper_bound_simple(5), "U_5 = 13")
    # representation used by bounds/known.py: R_n = (a - 7n sqrt z)/(18n)
    for n in range(4, 200):
        P = known.P(n)
        a = 36 * n * P + 25 * n * n - 43 * n - 12 * P + 36
        check(Fraction(a, 18 * n) == 2 * P + n - 2 - Fraction(2 * P, 3 * n) + Fraction(2, n) + Fraction(7 * (n - 1), 18),
              "R_n rational part", n)
    worst_slack = None
    for n in range(4, max_n + 1):
        U = known.strict_upper_bound(n)
        T2 = known.strict_upper_bound_simple(n)
        check(U <= T2, "(T) implies (T')", (n, U, T2))
        slack = T2 - U
        if worst_slack is None or slack < worst_slack[0]:
            worst_slack = (slack, n)
        if n % 2:      # odd: r_n = q(n-1) exactly -> z = 2(n-1)^2
            check(3 * n * n - 4 * n + 1 - 4 * known.P(n) == 2 * (n - 1) ** 2, "odd z", n)
        else:
            check(3 * n * n - 4 * n + 1 - 4 * known.P(n) == 2 * (n - 1) ** 2 - 1, "even z", n)
    check([known.strict_upper_bound(n) for n in (4, 7, 20, 100)] == [9, 27, 211, 5065], "control values U_n")
    return {"range": f"4<=n<={max_n}", "min_slack_T2_minus_U": worst_slack}


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cycle-max-n", type=int, default=9)
    ap.add_argument("--perm-max-n", type=int, default=7)
    ap.add_argument("--reflection-max-n", type=int, default=12)
    ap.add_argument("--scalar-max-n", type=int, default=100000)
    ap.add_argument("--seed", type=int, default=20260912)
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    if args.quick:
        args.cycle_max_n, args.perm_max_n, args.reflection_max_n, args.scalar_max_n = 7, 6, 8, 2000
    t0 = time.time()
    report = {"core_version": CORE_VERSION, "checker": "checks/check_C16.py",
              "source": "docs/incoming/LRX_STRICT_PROOF_RU-1.md", "args": vars(args)}
    timings = {}
    t = time.time(); report["cycles"] = run_cycles(args.cycle_max_n); timings["cycles"] = time.time() - t
    print("cycles:", report["cycles"], flush=True)
    t = time.time(); report["permutations_exhaustive"] = run_perms_exhaustive(args.perm_max_n); timings["perms"] = time.time() - t
    for n, s in report["permutations_exhaustive"].items():
        print(f"perms n={n}: {s['permutations']} perms, loss hist {s['loss_min_len_hist']}, "
              f"reduced {s['loss_min_reduced_hist']}, max (3.4)-B_n = {s['max_expr_3_4_minus_B_n']}", flush=True)
    t = time.time(); report["reflections_C24"] = run_reflections(args.reflection_max_n); timings["reflections"] = time.time() - t
    print("reflections OK up to n =", args.reflection_max_n, flush=True)
    spec = [(8, 30, None), (9, 20, None), (10, 20, None), (20, 10, None), (50, 4, None), (100, 3, None), (101, 3, None)]
    if args.quick:
        spec = [(8, 5, None), (20, 2, None)]
    t = time.time(); report["samples"] = run_samples(spec, args.seed); timings["samples"] = time.time() - t
    for n, s in report["samples"].items():
        worst = max(r["min_len"] - known.target_diameter(n) for r in s["rows"].values())
        print(f"samples n={n}: {s['cases']} cases, all words sort; max(min_c len - B_n) = {worst}", flush=True)
    bs = (1, 2) if args.quick else (1, 2, 3, 4)
    t = time.time(); report["family_C23"] = run_family_C23(bs); timings["family_C23"] = time.time() - t
    print("C23 family:", {n: (v["V"], v["T_min"], v["calE"]) for n, v in report["family_C23"].items()}, flush=True)
    t = time.time(); report["scalar"] = run_scalar(args.scalar_max_n); timings["scalar"] = time.time() - t
    print("scalar:", report["scalar"], flush=True)
    report["timings_seconds"] = {k: round(v, 2) for k, v in timings.items()}
    report["total_seconds"] = round(time.time() - t0, 2)
    report["result"] = "PASS"
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(os.path.join(OUT_DIR, "report.json"), "w") as f:
        json.dump(report, f, indent=1, default=str)
    print("PASS in", report["total_seconds"], "s")


if __name__ == "__main__":
    try:
        main()
    except CheckFailure as exc:
        print("FAIL:", exc)
        sys.exit(1)
