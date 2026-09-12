"""Construction N1 (docs/incoming/LRX_STRICT_PROOF_RU-1.md, sections 1-3) as a
reusable sorting algorithm.

CONSTRUCTION_VERSION = "strict_upper-1.0". Core: oracle-1.0 (oracle/moves.py).

The word for a permutation pi and a shift c is exactly the word of N1:
local words W_C of section 2 (carrier choice 2.1, contraction 2.2, blocks
(2.4), disjoint cancellations 2.4) applied at the first visit of each
carrier by the external route of section 3, whose length is H(c).

Two things here are NOT part of N1 and are marked as variants (AGENTS.md
rule 16): the "carrier" external route (shortest walk 0 -> carriers -> c,
length R_c <= H(c)) and free reduction of the assembled word. The lemmas of
N1 are not transferred to them; they are only measured (experiments/).

Shift rule for `sorting_word`: "min_len" takes the shift with the smallest
actual word length (all n shifts are built; ties -> smallest c). This is a
different algorithm from "min over c of the bound" only in name here, because
in N1 the actual length equals the bound 2F_c - S_c + H(c) exactly.

The letter-by-letter correctness of every produced word is asserted with
the reference moves; this module never reads distance tables.
"""

import os
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))

from moves import apply_word, freely_reduce, identity, delta, CORE_VERSION  # noqa: E402

CONSTRUCTION_VERSION = "strict_upper-1.0"


class ConstructionError(Exception):
    pass


def _req(cond, msg, info=None):
    if not cond:
        raise ConstructionError(f"{msg}: {info!r}")


# ------------------------------------------------------------------ N1 section 1

def signed_step(a, b, n):
    """d = b - a mod n with -n/2 < d <= n/2 (antipodal step positive, rule 16)."""
    d = (b - a) % n
    if 2 * d > n:
        d -= n
    return d


def arc_edges(a, d, n):
    """Indices v of edges {v, v+1} traversed by the arc from a with signed step d."""
    if d > 0:
        return [(a + i) % n for i in range(d)]
    return [(a - 1 - i) % n for i in range(-d)]


def cycles_of(f):
    """Nontrivial cycles of the map f (list), each as a list of positions."""
    n = len(f)
    seen = [False] * n
    out = []
    for i in range(n):
        if seen[i] or f[i] == i:
            seen[i] = True
            continue
        cyc, j = [], i
        while not seen[j]:
            seen[j] = True
            cyc.append(j)
            j = f[j]
        out.append(cyc)
    return out


def cycle_data(cycle, n):
    """Steps, loads, F, M, t, E, S of an oriented nontrivial cycle (N1 section 1)."""
    k = len(cycle)
    steps = [signed_step(cycle[j], cycle[(j + 1) % k], n) for j in range(k)]
    loads = [0] * n
    for j in range(k):
        for v in arc_edges(cycle[j], steps[j], n):
            loads[v] += 1
    F = sum(abs(d) for d in steps)
    M = max(loads)
    t = sum(1 for j in range(k) if steps[j - 1] < 0 < steps[j])
    E = k - 2 * t
    S = 2 * M + E - 1
    antipodal_transp = (k == 2 and n % 2 == 0 and steps[0] == n // 2 and steps[1] == n // 2)
    return {"k": k, "steps": steps, "loads": loads, "F": F, "M": M, "t": t,
            "E": E, "S": S, "antipodal_transposition": antipodal_transp}


# ------------------------------------------------------------------ N1 section 2

def choose_carrier(cycle, cd, n):
    """Section 2.1: index j0 in the cycle of the carrier b (edge {b,b+1} has load M).
    Unidirectional cycle: j0 = 0. Mixed: b = start of the first maximal run of load M."""
    steps, loads, M = cd["steps"], cd["loads"], cd["M"]
    if len({d > 0 for d in steps}) == 1:
        return 0
    starts = [v for v in range(n) if loads[v] == M and loads[v - 1] < M]
    b = starts[0]
    _req(b in cycle, "carrier not in cycle", (cycle, b))
    j0 = cycle.index(b)
    _req(steps[j0 - 1] < 0 < steps[j0], "(2.1) carrier must be a -->+ vertex", (cycle, b))
    return j0


def local_word(cycle, n):
    """Local word W_C of N1 section 2 in carrier coordinates (carrier at slot 0).

    Returns dict: carrier b, word (after the chosen cancellations, length 2F-S),
    raw word (2.4), the state in carrier coordinates that the word sorts
    (tuple q with q[v_j] = v_{j+1}, v = positions relative to b), and cd.
    """
    cd = cycle_data(cycle, n)
    j0 = choose_carrier(cycle, cd, n)
    k = cd["k"]
    cyc = cycle[j0:] + cycle[:j0]
    b = cyc[0]
    v = [(a - b) % n for a in cyc]
    steps = cd["steps"][j0:] + cd["steps"][:j0]
    e = []
    for j in range(k):
        passes = 0 in arc_edges(v[j], steps[j], n)          # arc crosses edge {0,1}
        e.append(steps[j] - (1 if steps[j] > 0 else -1) if passes else steps[j])
    blocks = [list("XL" * x) if x > 0 else list("RX" * (-x)) for x in e]
    tokens = []                                              # (letter, block, pos)
    for j in range(k):
        tokens.extend((ch, j, i) for i, ch in enumerate(blocks[j]))
        if j < k - 1:
            tokens.append(("X", -1, j + 1))                  # inserted X for vertex a_{j+1}
    raw = "".join(t[0] for t in tokens)
    kill = set()
    for j in range(1, k):
        ins = next(i for i, t in enumerate(tokens) if t[1] == -1 and t[2] == j)
        dout, din = steps[j], steps[j - 1]
        if dout > 0:
            pair = (ins, ins + 1)
        elif din < 0:
            pair = (ins - 1, ins)
        else:
            continue
        _req(not (set(pair) & kill), "cancellations must be disjoint", cycle)
        kill.update(pair)
    word = "".join(t[0] for i, t in enumerate(tokens) if i not in kill)
    _req(len(word) == 2 * cd["F"] - cd["S"], "(2.6) local word length", (cycle, word))
    state = list(range(n))
    for j in range(k):
        state[v[j]] = v[(j + 1) % k]
    state = tuple(state)
    _req(apply_word(state, word) == identity(n), "(2.5) local word must sort its cycle state", (cycle, word))
    return {"carrier": b, "v": v, "e": e, "word": word, "raw": raw, "state": state, "cd": cd}


# ------------------------------------------------------------------ N1 section 3

def H_of(c, n):
    """Sufficient external cost of N1 section 3."""
    c %= n
    return n if c == 0 else n + delta(0, c, n) - 2


def n1_route(c, n):
    """External route of N1 section 3 from 0 to c visiting every vertex: list of
    visited positions (starting with 0) and the letters (len = H(c))."""
    c %= n
    if c == 0:
        return [(i % n) for i in range(n + 1)], "L" * n
    d = delta(0, c, n)
    dirn = 1 if c <= n // 2 else -1
    pos, letters, cur = [0], [], 0
    plan = [(dirn, d - 1), (-dirn, d - 1), (-dirn, n - d)]
    for step, count in plan:
        for _ in range(count):
            cur = (cur + step) % n
            pos.append(cur)
            letters.append("L" if step > 0 else "R")
    _req(cur == c and len(letters) == H_of(c, n), "N1 route", (c, n))
    return pos, "".join(letters)


def carrier_route(required, c, n):
    """VARIANT (not N1): shortest walk on the cycle Z_n from 0 to c that visits every
    position in `required`. Returns (positions, letters). Brute force over the gap
    of the point set left uncovered (the walk covers a contiguous arc)."""
    c %= n
    pts = sorted(set(required) | {0, c})
    best = None
    m = len(pts)
    for g in range(m):
        # uncovered gap: strictly between pts[g] and pts[(g+1) % m]; covered arc from
        # lo = pts[(g+1) % m] forward to hi = pts[g]
        lo = pts[(g + 1) % m]
        hi = pts[g]
        L = (hi - lo) % n                    # arc length in edges (0 if single point)
        if m == 1:
            L = 0
        x0 = (0 - lo) % n                    # coordinate of 0 inside the arc
        xc = (c - lo) % n
        _req(x0 <= L and xc <= L, "arc must contain 0 and c", (pts, g))
        # option A: go down to lo first, then up to hi, then back to c
        costA = x0 + L + (L - xc)
        planA = [(-1, x0), (1, L), (-1, L - xc)]
        # option B: go up to hi first, then down to lo, then up to c
        costB = (L - x0) + L + xc
        planB = [(1, L - x0), (-1, L), (1, xc)]
        for cost, plan in ((costA, planA), (costB, planB)):
            if best is None or cost < best[0]:
                best = (cost, plan)
    # wrapping walks: one full turn in one direction, then on to c
    for dirn in (1, -1):
        extra = (dirn * c) % n
        cost = n + extra
        if cost < best[0]:
            best = (cost, [(dirn, n + extra)])
    pos, letters, cur = [0], [], 0
    for step, count in best[1]:
        for _ in range(count):
            cur = (cur + step) % n
            pos.append(cur)
            letters.append("L" if step > 0 else "R")
    _req(cur == c and set(pts) <= set(pos), "carrier route", (required, c, n))
    return pos, "".join(letters)


def shift_word(pi, c, route="n1"):
    """Full sorting word of N1 for permutation pi and shift c (route="n1"), or the
    variant with the carrier route (route="carrier"). Returns dict with the word,
    per-cycle data and rule-15 statistics; the word is asserted to sort pi."""
    n = len(pi)
    c %= n
    f = [(pi[i] + c) % n for i in range(n)]
    locals_ = [local_word(cyc, n) for cyc in cycles_of(f)]
    words = {}
    for lw in locals_:
        _req(lw["carrier"] not in words, "carriers must be distinct", (pi, c))
        words[lw["carrier"]] = lw["word"]
    if route == "n1":
        pos, letters = n1_route(c, n)
    elif route == "carrier":
        pos, letters = carrier_route(list(words), c, n)
    else:
        raise ValueError(route)
    parts, done = [], set()
    for idx, b in enumerate(pos):
        if b in words and b not in done:
            parts.append(words[b])
            done.add(b)
        if idx < len(letters):
            parts.append(letters[idx])
    _req(done == set(words), "every carrier must be visited", (pi, c))
    word = "".join(parts)
    _req(apply_word(tuple(pi), word) == identity(n), "full word must sort pi", (pi, c, word))
    F = sum(lw["cd"]["F"] for lw in locals_)
    S = sum(lw["cd"]["S"] for lw in locals_)
    stats = {
        "c": c, "route": route, "n_cycles": len(locals_),
        "F_c": F, "S_c": S,
        "sumM": sum(lw["cd"]["M"] for lw in locals_),
        "E_c": sum(lw["cd"]["E"] for lw in locals_),
        "theta_c": sum(1 for lw in locals_ if lw["cd"]["antipodal_transposition"]),
        "H": H_of(c, n), "route_len": len(letters),
        "carriers": [lw["carrier"] for lw in locals_],
        "local_len": sum(len(lw["word"]) for lw in locals_),
        "len": len(word), "N_X": word.count("X"),
    }
    stats["N_rot"] = stats["len"] - stats["N_X"]
    if route == "n1":
        _req(stats["len"] == 2 * F - S + H_of(c, n), "(3.3) length", (pi, c))
    return {"word": word, "locals": locals_, "stats": stats}


def sorting_word(state, rule="min_len", route="n1"):
    """Interface of PLAN 9.8: the sorting word for `state` with the shift rule.
    rule="min_len": all n shifts are built, the shortest word wins (ties: smallest c).
    Returns (word, chosen_c)."""
    n = len(state)
    _req(sorted(state) == list(range(n)) and n >= 2, "state must be a permutation", state)
    if rule != "min_len":
        raise ValueError(rule)
    best = None
    for c in range(n):
        sw = shift_word(state, c, route)
        if best is None or len(sw["word"]) < len(best[0]):
            best = (sw["word"], c)
    return best


def _selftest():
    import itertools
    import random
    # carrier_route against BFS on (position, visited-set) for small n
    def bfs_route(required, c, n):
        req = frozenset(set(required) | {0, c})
        start = (0, frozenset([0]) & req)
        dist = {start: 0}
        frontier = [start]
        while frontier:
            nxt = []
            for (p, vis) in frontier:
                if p == c and vis == req:
                    return dist[(p, vis)]
                for step in (1, -1):
                    q = (p + step) % n
                    v2 = vis | ({q} & req)
                    if (q, v2) not in dist:
                        dist[(q, v2)] = dist[(p, vis)] + 1
                        nxt.append((q, v2))
            frontier = nxt
        raise AssertionError
    cnt = 0
    for n in range(2, 9):
        for c in range(n):
            for r in range(n + 1):
                for req in itertools.combinations(range(n), r):
                    pos, letters = carrier_route(list(req), c, n)
                    assert len(letters) == bfs_route(req, c, n), (req, c, n)
                    assert len(letters) <= H_of(c, n)
                    cnt += 1
    print("carrier_route == BFS on", cnt, "cases, 2<=n<=8")
    rng = random.Random(1)
    for n in range(4, 9):
        for _ in range(50):
            p = list(range(n)); rng.shuffle(p)
            for c in range(n):
                shift_word(tuple(p), c, "n1"); shift_word(tuple(p), c, "carrier")
    print("shift_word OK on random inputs 4<=n<=8;", CORE_VERSION, CONSTRUCTION_VERSION)


if __name__ == "__main__":
    _selftest()
