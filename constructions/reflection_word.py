"""Two-arc reversal: an explicit word sorting every reflection pi_h(i) = h - i mod n
in B_n - delta_n(h, 1) moves (C32, docs/proofs/C32_two_arc_reversal.md).

Circle model (PROBLEM section 2): circle[pos] = element, head at position 0
initially; a shift c means that element e must end at circle position e + c and
the head must end at c.  For pi_h the required position map is
    pos i  ->  h' - i  (mod n),   h' = h + c,
a reflection of Z_n with axis through the points h'/2 and h'/2 + n/2 (vertices
when the corresponding number is an integer, edge midpoints otherwise).  The
two arcs centred at the two axis points are each mapped onto themselves
reversed, so the reflection is the product of two LINEAR reversals of arcs of
lengths a + b = n (a, b odd for a vertex-centred arc, even for an edge-centred
arc).  A linear reversal of an arc of length a needs a(a-1)/2 adjacent swaps,
and the "zigzag" schedules below make consecutive swaps at adjacent edges, so
the head pays one rotation between consecutive swaps.

Edges of the circle: e = {e, e+1}; the head sees edge e when head == e.

Zigzag schedules for a linear reversal of the arc of positions lo..lo+a-1
(returned as a list of edges, consecutive entries always adjacent):
  * "out": centre-out, starts at the central edge (even a) or at the edge
    right of the central vertex (odd a), ends at the topmost edge lo+a-2;
  * "out_mirror": mirror image, ends at the lowest edge lo;
  * "in": reverse of "out_mirror": starts at the lowest edge, ends at the
    centre (edge lo+a/2-1 for even a, edge lo+(a-1)/2-1 for odd a);
  * "in_mirror": reverse of "out": starts at the topmost edge, ends at the
    centre (edge lo+a/2-1 for even a, edge lo+(a-1)/2 for odd a).
Any of them reverses the arc (a reversal is an involution, so the reversed
edge sequence is again a reversal; the mirror image reverses the arc by
symmetry).

reflection_word(n, h) chooses h', c, the arc lengths and the schedules as in the
proof (by n mod 4 and the position of h), builds the word, checks it with the
reference moves and returns rule-15 statistics.  The word has
    N_X = floor((n-1)^2/4),  N_rot = floor(n^2/4) - delta_n(h, 1),
    len = B_n - delta_n(h, 1).
The function never reads distance tables.

policy_like(n, hp, first) builds the variant that the "horizon" policy
(candidate_sweep sweep-1.2, c = 0, dir = +1) produces on reflections: both
arcs centre-out, the head walking from the end of the first arc to the centre
of the second; used by checks/check_C32.py to identify the policy's words.
"""

import os
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))

from moves import apply_word, identity, delta, CORE_VERSION  # noqa: E402

CONSTRUCTION_VERSION = "two_arc-1.0"


def zigzag_out(lo, a):
    """Centre-out linear reversal of positions lo..lo+a-1 (a >= 1): list of edges.

    Even a = 2s: round 0 swaps the central edge m = lo+s-1; round t = 1..s-1
    walks down from m+t to m-t (the top element of the block moves to its
    bottom) and up from m-t+1 to m+t (the bottom element moves to the top).
    Odd a = 2s+1, centre vertex v = lo+s: round t = 0..s-1 walks down from v+t
    to v-t-1 and up from v-t to v+t.  Consecutive edges differ by exactly 1."""
    seq = []
    if a % 2 == 0:
        s = a // 2
        m = lo + s - 1
        seq.append(m)
        for t in range(1, s):
            seq.extend(range(m + t, m - t - 1, -1))
            seq.extend(range(m - t + 1, m + t + 1))
    else:
        s = a // 2
        v = lo + s
        for t in range(0, s):
            seq.extend(range(v + t, v - t - 2, -1))
            seq.extend(range(v - t, v + t + 1))
    return seq


def zigzag(lo, a, mode):
    """Edge sequence for the linear reversal of positions lo..lo+a-1.
    mode in {"out", "out_mirror", "in", "in_mirror"} (see module docstring)."""
    if mode == "out":
        return zigzag_out(lo, a)
    if mode == "in_mirror":
        return list(reversed(zigzag_out(lo, a)))
    # mirror: position p -> 2*lo + a - 1 - p, edge e -> 2*lo + a - 2 - e
    mirrored = [2 * lo + a - 2 - e for e in zigzag_out(lo, a)]
    if mode == "out_mirror":
        return mirrored
    if mode == "in":
        return list(reversed(mirrored))
    raise ValueError(mode)


def edges_to_word(edges, n, start=0, end=None):
    """Word that moves the head from `start` along the shorter arc to each edge in
    turn, swapping there, and finally to `end` (if given) by the shorter arc."""
    word = []
    head = start
    for e in edges:
        e %= n
        fwd = (e - head) % n
        if fwd <= n - fwd:
            word.append("L" * fwd)
        else:
            word.append("R" * (n - fwd))
        word.append("X")
        head = e
    if end is not None:
        end %= n
        fwd = (end - head) % n
        word.append("L" * fwd if fwd <= n - fwd else "R" * (n - fwd))
    return "".join(word)


def apply_edges(circle, edges):
    n = len(circle)
    c = list(circle)
    for e in edges:
        a, b = e % n, (e + 1) % n
        c[a], c[b] = c[b], c[a]
    return c


def two_arc_edges(n, twoA, mode1, mode2, a=None):
    """Edge sequence reversing the circle about the axis through the point with
    doubled coordinate twoA (a vertex twoA/2 if twoA is even, the midpoint of the
    edge {(twoA-1)/2, (twoA+1)/2} if odd) and its antipode twoA + n.  The arc
    centred at twoA/2 is reversed first with schedule mode1, the arc centred at
    the antipode with schedule mode2.  The axis parameter is h' = twoA mod n
    (position map i -> h' - i).  `a` is the length of the first arc (default:
    the balanced choice of the proof; a vertex-centred arc is odd, an
    edge-centred arc even).  Returns (edges, lo1, a, lo2, b)."""
    twoA %= 2 * n
    twoB = (twoA + n) % (2 * n)
    if a is None:
        if n % 2 == 0:
            a = n // 2
            if (a % 2 == 1) != (twoA % 2 == 0):
                a -= 1
        else:
            m = n // 2
            odd_len = m if m % 2 == 1 else m + 1
            a = odd_len if twoA % 2 == 0 else n - odd_len
    b = n - a
    if (a % 2 == 1) != (twoA % 2 == 0) or (b % 2 == 1) != (twoB % 2 == 0):
        raise ValueError(f"parity of arcs: n={n} twoA={twoA} a={a} b={b}")
    lo1 = (twoA - (a - 1)) // 2          # first position of the first arc
    lo2 = (twoB - (b - 1)) // 2
    # the two arcs are complementary on the circle
    assert (lo1 + a) % n == lo2 % n and (lo2 + b) % n == lo1 % n
    edges = [e % n for e in zigzag(lo1, a, mode1)] + [e % n for e in zigzag(lo2, b, mode2)]
    return edges, lo1 % n, a, lo2 % n, b


def _stats(word, pi, n, c, extra):
    if apply_word(tuple(pi), word) != identity(n):
        raise AssertionError(f"word does not sort: n={n} pi={pi} c={c}")
    nx = word.count("X")
    out = {"word": word, "len": len(word), "N_X": nx, "N_rot": len(word) - nx, "c": c}
    out.update(extra)
    return out


def reflection_word(n, h):
    """Word of length B_n - delta_n(h, 1) sorting pi_h(i) = h - i mod n (n >= 4).
    Case analysis of the proof (docs/proofs/C32_two_arc_reversal.md, section 4)."""
    h %= n
    pi = tuple((h - i) % n for i in range(n))
    if n % 2 == 0:
        # axis point at the head: edge {0, 1} (h' = 1) when n = 0 mod 4, vertex 0
        # (h' = 0) when n = 2 mod 4; first arc centre-out, second boundary-in.
        twoA = 1 if n % 4 == 0 else 0
        c = (twoA - h) % n
        edges, lo1, a, lo2, b = two_arc_edges(n, twoA, "out", "in")
        word = edges_to_word(edges, n, 0, c)
        return _stats(word, pi, n, c, {"h": h, "h_prime": twoA % n, "arcs": (a, b), "case": "even"})
    m = n // 2
    if 1 <= h <= m:
        # case A: axis h' = 0; vertex arc around 0 first, centre-out (starts at
        # head 0); edge arc around {m, m+1} boundary-in, ends at head m.
        c = (-h) % n
        edges, lo1, a, lo2, b = two_arc_edges(n, 0, "out", "in")
        word = edges_to_word(edges, n, 0, c)
        return _stats(word, pi, n, c, {"h": h, "h_prime": 0, "arcs": (a, b), "case": "odd_A"})
    if h == 0 or m + 2 <= h <= n - 1:
        # case B: axis through the vertex v = h - 1, h' = 2v; the edge arc around
        # {v+m, v+m+1} first, centre-out (starts at head v+m); the vertex arc
        # boundary-in, ending at edge {v-1, v}, head v-1 = c.
        v = (h - 1) % n
        c = (2 * v - h) % n
        assert c == (v - 1) % n
        edges, lo1, a, lo2, b = two_arc_edges(n, 2 * v + n, "out", "in")
        word = edges_to_word(edges, n, 0, c)
        return _stats(word, pi, n, c, {"h": h, "h_prime": (2 * v) % n, "arcs": (a, b), "case": "odd_B"})
    # h = m + 1: mirror of h = m + 2 (conjugation by sigma_n swaps L and R and
    # maps pi_h to pi_{2-h}, PROBLEM section 4.3).
    r = reflection_word(n, (2 - h) % n)
    word = r["word"].translate(str.maketrans("LR", "RL"))
    return _stats(word, pi, n, (-r["c"]) % n,
                  {"h": h, "h_prime": (2 - r["h_prime"]) % n, "arcs": r["arcs"], "case": "odd_mirror"})


def policy_like(n, twoA):
    """Variant reproduced by the horizon policy on reflections (c = 0, dir = +1):
    both arcs centre-out, first the arc centred at doubled coordinate twoA."""
    hp = twoA % n
    pi = tuple((hp - i) % n for i in range(n))
    edges, lo1, a, lo2, b = two_arc_edges(n, twoA, "out", "out")
    word = edges_to_word(edges, n, 0, 0)
    return _stats(word, pi, n, 0, {"h": hp, "h_prime": hp, "arcs": (a, b), "edges": edges})


def _selftest():
    from moves import sigma, rev
    for n in range(4, 41):
        B = n * (n - 1) // 2
        for h in range(n):
            r = reflection_word(n, h)
            assert r["len"] == B - delta(h, 1, n), (n, h, r["len"])
            assert r["N_X"] == (n - 1) ** 2 // 4, (n, h, r["N_X"])
            assert r["N_rot"] == n * n // 4 - delta(h, 1, n), (n, h, r["N_rot"])
        assert reflection_word(n, 1)["len"] == B
        assert reflection_word(n, n - 1)["len"] == B - 2
        assert apply_word(sigma(n), reflection_word(n, 1)["word"]) == identity(n)
        assert apply_word(rev(n), reflection_word(n, n - 1)["word"]) == identity(n)
    print("reflection_word selftest ok", CONSTRUCTION_VERSION, CORE_VERSION)


if __name__ == "__main__":
    _selftest()
