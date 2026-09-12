"""N1 construction (docs/incoming/LRX_STRICT_PROOF_RU-1.md, §§1-3), written
from the proof text as an independent implementation (PLAN.md §3.1 item 7,
§9.10a). It is NOT the author's ``verify_lrx.py`` (not supplied) and it does
not import ``docs/incoming/lrx_2n3_check.py``.

Conventions are those of PROBLEM.md §3 (0-indexed, words left to right,
moves from oracle/moves.py). For a permutation pi and a shift c the word for
f_c(i) = pi(i) + c (mod n) is assembled exactly as in N1 §3: local words
(N1 §2) at the carriers, the carriers visited by the route of N1 §3 from
head position 0 to head position c.

Statistics follow AGENTS.md rule 15: for every (pi, c) we record F_c, S_c,
sum M(C), E_c, theta_c, H(c), the actual route length, local word lengths,
the full word length before and after free reduction, N_X and N_rot.
"""

import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "oracle"))
from moves import apply_word, freely_reduce, delta, word_stats  # noqa: E402

CONSTRUCTION_VERSION = "strict_upper-1.0 (N1 sections 1-3, oracle-1.0)"


# ---------------------------------------------------------------- N1 §1 ---

def signed_step(a, b, n):
    """Signed shortest step a -> b on Z_n with -n/2 < d <= n/2 (N1 §1,
    AGENTS.md rule 16: the antipodal step is positive)."""
    d = (b - a) % n
    return d if 2 * d <= n else d - n


def arc_edges(a, d, n):
    """Edges {v, v+1} (encoded by v) used by the signed step d from a."""
    if d > 0:
        return [(a + j) % n for j in range(d)]
    return [(a - 1 - j) % n for j in range(-d)]


def cycle_stats(cycle, n):
    """F, M, t, E, S and the edge loads of one nontrivial cycle (N1 §1)."""
    k = len(cycle)
    assert k >= 2
    ds = [signed_step(cycle[j], cycle[(j + 1) % k], n) for j in range(k)]
    loads = [0] * n
    for a, d in zip(cycle, ds):
        for e in arc_edges(a, d, n):
            loads[e] += 1
    F = sum(abs(d) for d in ds)
    M = max(loads)
    t = sum(1 for j in range(k) if ds[j - 1] < 0 < ds[j])
    t_minus = sum(1 for j in range(k) if ds[j - 1] > 0 > ds[j])
    E = k - 2 * t
    S = 2 * M + E - 1
    antipodal_transposition = (k == 2 and n % 2 == 0
                               and ds[0] == ds[1] == n // 2)
    return {
        "k": k, "ds": ds, "loads": loads, "F": F, "M": M, "t": t,
        "t_minus": t_minus, "E": E, "S": S,
        "mixed": any(d > 0 for d in ds) and any(d < 0 for d in ds),
        "antipodal_transposition": antipodal_transposition,
    }


def positive_support_is_one_arc(loads):
    """True iff {v : loads[v] > 0} is a single circular arc (or everything).

    N1 §1 (last paragraph) and §7.1 use that the positive support is
    connected."""
    n = len(loads)
    pos = [v for v in range(n) if loads[v] > 0]
    if len(pos) in (0, n):
        return True
    # count maximal runs of positive load on the circle
    runs = sum(1 for v in range(n) if loads[v] > 0 and loads[v - 1] == 0)
    return runs == 1


# ---------------------------------------------------------------- N1 §2 ---

def choose_carrier(cycle, st, n):
    """Index j into `cycle` of the carrier b = a_j (N1 §2.1).

    Mixed case: start of a maximal plateau of load M; by (1.1) it is a
    vertex of the cycle with a (- -> +) turn. One-sign case: a_0."""
    if not st["mixed"]:
        return 0
    M, loads = st["M"], st["loads"]
    cands = [j for j, a in enumerate(cycle)
             if loads[a] == M and loads[(a - 1) % n] < M]
    assert cands, "no maximal-load plateau start on the cycle (N1 §2.1)"
    j = cands[0]
    ds = st["ds"]
    assert ds[j - 1] < 0 < ds[j], "plateau start is not a (- -> +) turn"
    return j


def local_word(cycle, n):
    """Local word W_C of N1 §2 for one nontrivial cycle.

    Returns dict with: carrier b, coordinates v_j, contracted steps e_j,
    raw word (before cancellations), word after the disjoint XX
    cancellations of §2.4 (still a literal subsequence-deletion of the raw
    word), the number of cancellations, and the freely reduced word.
    All words are in the carrier's coordinate system (head at position b).
    """
    st = cycle_stats(cycle, n)
    j0 = choose_carrier(cycle, st, n)
    rot = list(cycle[j0:]) + list(cycle[:j0])
    k = len(rot)
    b = rot[0]
    ds = [signed_step(rot[j], rot[(j + 1) % k], n) for j in range(k)]
    if st["mixed"]:
        assert ds[0] > 0 and ds[-1] < 0, "(2.1) violated"
    # contraction of the edge {b, b+1} (edge index b): (2.2)
    es = []
    for j in range(k):
        uses = b in arc_edges(rot[j], ds[j], n)
        es.append(ds[j] - (1 if ds[j] > 0 else -1) if uses else ds[j])
    assert sum(abs(e) for e in es) == st["F"] - st["M"], "(2.2) violated"
    assert sum(es) % (n - 1) == 0, "(2.3) violated"
    for j in range(1, k - 1):
        assert es[j] != 0 and (es[j] > 0) == (ds[j] > 0), "§2.2 inner step"
    # blocks: U^e = (XL)^e for e > 0, (RX)^{-e} for e < 0 ; X between (2.4)
    blocks = []
    for e in es:
        blocks.append(["X", "L"] * e if e > 0 else ["R", "X"] * (-e))
    raw = []
    for j, blk in enumerate(blocks):
        raw.extend(blk)
        if j < k - 1:
            raw.append("X")
    raw = "".join(raw)
    assert len(raw) == 2 * (st["F"] - st["M"]) + k - 1

    # disjoint cancellations of §2.4, applied literally as deletions
    # position bookkeeping: block j occupies [start[j], start[j]+len)
    # followed by inserted X_j (for j < k-1)
    starts, pos = [], 0
    for j, blk in enumerate(blocks):
        starts.append(pos)
        pos += len(blk) + (1 if j < k - 1 else 0)
    delete = set()
    cancels = 0
    for j in range(1, k):            # inserted X_j sits after block j-1
        xj = starts[j - 1] + len(blocks[j - 1])
        if ds[j] > 0:                 # outgoing positive: X_j with first X of block j
            assert blocks[j] and blocks[j][0] == "X"
            delete.update({xj, starts[j]})
            cancels += 1
        elif ds[j - 1] < 0:           # incoming negative: last X of block j-1 with X_j
            assert blocks[j - 1] and blocks[j - 1][-1] == "X"
            delete.update({xj - 1, xj})
            cancels += 1
        # else (+ -> -): no cancellation
    assert len(delete) == 2 * cancels, "cancellation pairs overlap"
    assert cancels == k - 1 - st["t"], "number of cancellations != k-1-t"
    word = "".join(ch for i, ch in enumerate(raw) if i not in delete)
    assert len(word) == 2 * st["F"] - st["S"], "(2.6) not attained exactly"
    reduced = freely_reduce(word)

    vs = [(a - b) % n for a in rot]
    return {
        "stats": st, "carrier": b, "rot": rot, "vs": vs, "es": es,
        "raw": raw, "word": word, "reduced": reduced, "cancels": cancels,
    }


def local_word_action_ok(lw, n):
    """Check (2.5): the position permutation of the local word equals
    C'^{-1} in the carrier's coordinates (C' = (0 v_1 ... v_{k-1}))."""
    q = apply_word(tuple(range(n)), lw["word"])     # q[i] = g(i)
    vs = lw["vs"]
    k = len(vs)
    expect = list(range(n))
    for j in range(k):
        expect[vs[(j + 1) % k]] = vs[j]              # C'^{-1}(v_{j+1}) = v_j
    return list(q) == expect


# ---------------------------------------------------------------- N1 §3 ---

def H(c, n):
    c %= n
    return n if c == 0 else n + delta(0, c, n) - 2


def route(n, c):
    """External route of N1 §3 from head position 0 to c visiting every
    vertex. Returns a list of (move, head_position_after_move)."""
    c %= n
    out, pos = [], 0
    if c == 0:
        for _ in range(n):
            pos = (pos + 1) % n
            out.append(("L", pos))
        return out
    d = signed_step(0, c, n)
    fwd, back = ("L", "R") if d > 0 else ("R", "L")
    step = 1 if d > 0 else -1
    for _ in range(abs(d) - 1):
        pos = (pos + step) % n
        out.append((fwd, pos))
    for _ in range(abs(d) - 1):
        pos = (pos - step) % n
        out.append((back, pos))
    for _ in range(n - abs(d)):
        pos = (pos - step) % n
        out.append((back, pos))
    assert pos == c and len(out) == H(c, n)
    assert len({p for _, p in out} | {0}) == n
    return out


def nontrivial_cycles(f):
    n = len(f)
    seen = [False] * n
    cycles = []
    for s in range(n):
        if seen[s]:
            continue
        cyc, x = [], s
        while not seen[x]:
            seen[x] = True
            cyc.append(x)
            x = f[x]
        if len(cyc) > 1:
            cycles.append(cyc)
    return cycles


def shift_word(pi, c):
    """Full word of N1 §3 for permutation pi and shift c, with statistics."""
    n = len(pi)
    f = [(pi[i] + c) % n for i in range(n)]
    cycles = nontrivial_cycles(f)
    locals_ = [local_word(C, n) for C in cycles]
    by_carrier = {lw["carrier"]: lw for lw in locals_}
    assert len(by_carrier) == len(locals_)
    parts = []
    route_moves = 0
    used = set()
    if 0 in by_carrier:                 # head starts at 0
        parts.append(by_carrier[0]["word"])
        used.add(0)
    for mv, pos in route(n, c):
        parts.append(mv)
        route_moves += 1
        if pos in by_carrier and pos not in used:   # first visit only
            parts.append(by_carrier[pos]["word"])
            used.add(pos)
    assert used == set(by_carrier), "route missed a carrier"
    word = "".join(parts)
    st = [lw["stats"] for lw in locals_]
    F_c = sum(s["F"] for s in st)
    S_c = sum(s["S"] for s in st)
    E_c = sum(s["E"] for s in st)
    sumM = sum(s["M"] for s in st)
    theta = sum(1 for s in st if s["antipodal_transposition"])
    local_len = sum(len(lw["word"]) for lw in locals_)
    reduced = freely_reduce(word)
    length, n_x, n_rot = word_stats(word)
    rlength, rn_x, rn_rot = word_stats(reduced)
    return {
        "c": c, "word": word, "reduced": reduced,
        "cycles": cycles, "locals": locals_,
        "F_c": F_c, "S_c": S_c, "E_c": E_c, "sumM": sumM, "theta_c": theta,
        "H": H(c, n), "route_len": route_moves,
        "local_len": local_len, "local_bound": 2 * F_c - S_c,
        "bound": 2 * F_c - S_c + H(c, n),
        "len": length, "N_X": n_x, "N_rot": n_rot,
        "reduced_len": rlength, "reduced_N_X": rn_x, "reduced_N_rot": rn_rot,
    }


def all_shifts(pi):
    return [shift_word(pi, c) for c in range(len(pi))]


# ------------------------------------------------- N1 §§4-6 statistics ---

def window(i, n):
    return {(i + j) % n for j in range(1, n // 2 + 1)}


def reflection_errors(pi):
    """e_i, V, x of N1 §4 (tau = -pi)."""
    n = len(pi)
    tau = [(-v) % n for v in pi]
    es = []
    for i in range(n):
        img = {tau[j] for j in window(i, n)}
        es.append(len(img ^ window(tau[i], n)))
    V = sum(es)
    return {"tau": tau, "e": es, "V": V, "x": Fraction(V, n)}


def T(pi, h):
    n = len(pi)
    return sum(delta((pi[i] + i) % n, h, n) for i in range(n))


def antipodal_pairs_preserved(pi):
    """B of N1 §4: unordered pairs antipodal before and after pi (even n)."""
    n = len(pi)
    if n % 2:
        return 0
    return sum(1 for i in range(n // 2)
               if (pi[i] - pi[i + n // 2]) % n == n // 2)


def perm_summary(pi, shifts=None):
    """Averaged N1 quantities for one permutation (exact fractions)."""
    n = len(pi)
    if shifts is None:
        shifts = all_shifts(pi)
    P = n * n // 4
    avg = lambda key: Fraction(sum(s[key] for s in shifts), n)  # noqa: E731
    re = reflection_errors(pi)
    Ts = [T(pi, h) for h in range(n)]
    return {
        "n": n, "P": P,
        "F_bar": avg("F_c"), "H_bar": avg("H"), "s": avg("S_c"),
        "mu": avg("sumM"), "E_bar": avg("E_c"), "theta_bar": avg("theta_c"),
        "calE": sum(s["E_c"] for s in shifts),
        "K_bar": Fraction(sum(len(s["cycles"]) for s in shifts), n),
        "V": re["V"], "x": re["x"], "e": re["e"],
        "B": antipodal_pairs_preserved(pi),
        "T": Ts, "T_min": min(Ts),
        "avg_bound": avg("bound"),
        "min_bound": min(s["bound"] for s in shifts),
        "min_len": min(s["len"] for s in shifts),
        "min_reduced_len": min(s["reduced_len"] for s in shifts),
    }


if __name__ == "__main__":
    # PROBLEM.md §3 mandatory convention example and a small demonstration
    assert apply_word((2, 0, 1), "XL") == (2, 1, 0)
    from moves import sigma  # noqa: E402
    for n in (4, 5, 6, 7):
        pi = sigma(n)
        sh = all_shifts(pi)
        for s in sh:
            assert apply_word(pi, s["word"]) == tuple(range(n))
        best = min(sh, key=lambda s: (s["len"], s["c"]))
        print(f"n={n} sigma_n: min_c len={best['len']} (c={best['c']}, "
              f"N_X={best['N_X']}, N_rot={best['N_rot']}), "
              f"min_c bound={min(s['bound'] for s in sh)}, "
              f"B_n={n*(n-1)//2}")
