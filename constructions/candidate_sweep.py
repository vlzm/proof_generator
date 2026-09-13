"""Candidate construction "sweep": joint processing of all cycles by head sweeps.

Motivation (session 5, data/runs/geodesic_structure/): on the affine inputs where
N1 loses Theta(n), every geodesic walks the head once forward around the circle
and once backward, pairing forward-moving and backward-moving elements of
DIFFERENT cycles of f_c in one swap.  N1 instead sweeps each cycle separately
and pays the route between carriers.

Algorithm (CANDIDATE_VERSION "sweep-1.0"), for a permutation pi, a shift c and
a first direction dir in {+1, -1}:

  * circle model (PROBLEM section 2): circle[pos] = element, head at position
    0; element e must end at circle position (e + c) mod n; rem(pos) = signed
    shortest step from pos to the target of circle[pos] (antipodal step
    positive, rule 16).
  * a pass walks the head in direction dir (L for +1, R for -1).  At the edge
    {h, h+1} under the head, with x = circle[h], y = circle[h+1], the swap X
    is made iff
        gain = [rem x > 0] + [rem y < 0] - [rem x < 0] - [rem y > 0] >= thr
    (thr = 1: also carry an element past a bystander that does not mind, i.e.
    a "single" swap; thr = 2: only double swaps).  The pass ends when no
    element strictly ahead of the head (within the n-1 remaining steps) is
    displaced, or when everything is in place.
  * passes alternate direction until Phi = sum |rem| = 0; then the head goes
    to c by the shorter arc (the word must end with the head at c so that the
    array reads as the identity).
  * safety: at most 4n passes, otherwise ConstructionError for that (c, dir).

Mode "nearest" (sweep-1.1): instead of full passes, after each swap (or at an
unproductive edge) the head walks towards the nearest edge where a swap with
gain >= thr is possible in the current configuration (ties keep the current
direction).  On sigma_n this reproduces the short oscillations of the
geodesics; mode "pass" reproduces the two long sweeps of the geodesics on
affine inputs.  best_sweep takes the shortest word over both modes.

Mode "horizon" (sweep-1.2, H12 candidate (a)): at every decision point the
head chooses between the two directions by a bounded look-ahead.  A branch
walks to the nearest productive edge (gain >= thr) in its direction, swaps
there, and recurses; after `horizon` swaps (or Phi = 0) the leaf is scored
    f = rotations + swaps so far + weight * Phi_leaf
(Phi_leaf/2 is a lower bound on the remaining swaps; weight 1 charges one
rotation per remaining swap as well).  The first step of the best branch is
executed; ties keep the current direction.  horizon = 1 with weight 0 is
mode "nearest".

Deadlock (sweep-1.3, `fallback=True`): under a fixed shift c a cycle whose
elements all travel in the same direction along their shortest arcs has no
edge with gain >= 1 (e.g. n = 7, pi = (0,1,3,6,4,5,2), c = 0, cycle (2 3 6):
steps +1, +3, +3).  Without the fallback the (c, dir, thr) attempt raises
ConstructionError and best_sweep relies on another shift (at n = 6, 7 more
than half of all (pi, c) pairs are dead).  With fallback=True, when no edge
with gain >= thr exists in either direction, the displaced element with the
largest |rem| is marked as travelling the long way round (its direction of
travel is fixed to the opposite sign, a persistent flag: the role of N1's
carrier; extra distance n - 2|rem|), and
the look-ahead continues with the modified rem.  Gain-0 swaps are never
used (they oscillate).  fallback=False reproduces sweep-1.2 exactly.

`sweep_word(pi, c, dir, thr)` returns the word and rule-15 statistics and
asserts with the reference moves that the word sorts pi.  `best_sweep(pi)`
minimises the (unreduced) length over c, dir and thr in {1, 2}.
This module never reads distance tables.
"""

import os
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))

from moves import apply_word, freely_reduce, identity, delta, CORE_VERSION  # noqa: E402

CANDIDATE_VERSION = "sweep-1.3"


class ConstructionError(Exception):
    pass


def signed_step(a, b, n):
    d = (b - a) % n
    if 2 * d > n:
        d -= n
    return d


def _productive_ahead(circle, target, head, d, n, thr, rem):
    """Distance (in head steps, 1..n-1) to the nearest edge in direction d where a
    swap with gain >= thr is possible in the current configuration; None if none."""
    for k in range(1, n):
        h = (head + d * k) % n
        a, b = h, (h + 1) % n
        ra, rb = rem(a), rem(b)
        if (ra > 0) + (rb < 0) - (ra < 0) - (rb > 0) >= thr:
            return k
    return None


def sweep_word(pi, c, dirn=1, thr=1, mode="pass", max_passes=None, horizon=3, weight=1.0,
               fallback=False):
    """mode="pass": alternating full passes (sweep-1.0).  mode="nearest": after each
    swap (or when the current edge is unproductive) the head goes towards the
    nearest productive edge, ties keep the current direction (sweep-1.1)."""
    n = len(pi)
    c %= n
    if max_passes is None:
        max_passes = 4 * n
    circle = list(pi)
    target = [0] * n
    for e in range(n):
        target[e] = (e + c) % n
    head = 0
    word = []
    phi = sum(delta(i, target[circle[i]], n) for i in range(n))

    def rem(pos):
        return signed_step(pos, target[circle[pos]], n)

    passes = 0
    flips = 0
    if mode == "pass":
        while phi > 0:
            if passes >= max_passes:
                raise ConstructionError(f"too many passes: {pi} c={c} dir={dirn} thr={thr}")
            d = dirn if passes % 2 == 0 else -dirn
            for step in range(n - 1):
                a, b = head, (head + 1) % n
                ra, rb = rem(a), rem(b)
                gain = (ra > 0) + (rb < 0) - (ra < 0) - (rb > 0)
                if gain >= thr:
                    circle[a], circle[b] = circle[b], circle[a]
                    word.append("X")
                    phi = sum(delta(i, target[circle[i]], n) for i in range(n))
                    if phi == 0:
                        break
                # anything displaced strictly ahead within the remaining steps?
                ahead = False
                for k in range(1, n - step):
                    p = (head + d * k) % n
                    if rem(p) != 0:
                        ahead = True
                        break
                if not ahead:
                    break
                head = (head + d) % n
                word.append("L" if d > 0 else "R")
            passes += 1
    elif mode == "nearest":
        d = dirn
        guard = 0
        while phi > 0:
            guard += 1
            if guard > 8 * n * n:
                raise ConstructionError(f"nearest mode does not terminate: {pi} c={c}")
            a, b = head, (head + 1) % n
            ra, rb = rem(a), rem(b)
            gain = (ra > 0) + (rb < 0) - (ra < 0) - (rb > 0)
            if gain >= thr:
                circle[a], circle[b] = circle[b], circle[a]
                word.append("X")
                phi = sum(delta(i, target[circle[i]], n) for i in range(n))
                continue
            kf = _productive_ahead(circle, target, head, 1, n, thr, rem)
            kb = _productive_ahead(circle, target, head, -1, n, thr, rem)
            if kf is None and kb is None:
                raise ConstructionError(f"no productive edge: {pi} c={c} thr={thr}")
            if kb is None or (kf is not None and (kf < kb or (kf == kb and d > 0))):
                nd = 1
            else:
                nd = -1
            if nd != d:
                passes += 1
            d = nd
            head = (head + d) % n
            word.append("L" if d > 0 else "R")
    elif mode == "horizon":
        d = dirn
        guard = 0
        H = max(1, int(horizon))
        flips = 0

        long_arc = {}   # element -> forced direction (+1/-1) of travel (fallback)

        def rem_of(circ, pos):
            e = circ[pos]
            sgn = long_arc.get(e)
            if sgn is None:
                return signed_step(pos, target[e], n)
            if sgn > 0:
                return (target[e] - pos) % n
            return -((pos - target[e]) % n)

        def phi_of(circ):
            return sum(abs(rem_of(circ, i)) for i in range(n))

        def nearest_edge(circ, hd, dd):
            for k in range(0, n):
                h = (hd + dd * k) % n
                ra, rb = rem_of(circ, h), rem_of(circ, (h + 1) % n)
                if (ra > 0) + (rb < 0) - (ra < 0) - (rb > 0) >= thr:
                    return k
            return None

        def search(circ, hd, ph, depth, cost, cur_d):
            """Best f over branches of remaining depth; returns (f, first_step)."""
            if ph == 0 or depth == 0:
                return cost + weight * ph, None
            best = None
            for dd in (cur_d, -cur_d):
                k = nearest_edge(circ, hd, dd)
                if k is None:
                    continue
                h = (hd + dd * k) % n
                c2 = list(circ)
                a, b = h, (h + 1) % n
                c2[a], c2[b] = c2[b], c2[a]
                ph2 = phi_of(c2)
                nd = dd if k > 0 else cur_d
                f, _ = search(c2, h, ph2, depth - 1, cost + k + 1, nd)
                if best is None or f < best[0]:
                    best = (f, dd, k)
            if best is None:
                return cost + weight * ph, None
            return best[0], (best[1], best[2])

        while phi > 0:
            guard += 1
            if guard > 8 * n * n:
                raise ConstructionError(f"horizon mode does not terminate: {pi} c={c}")
            f, step = search(circle, head, phi, H, 0, d)
            if step is None:
                if not fallback:
                    raise ConstructionError(f"no productive edge: {pi} c={c} thr={thr}")
                # deadlock: route the most displaced element the long way round
                cand = [(abs(rem_of(circle, i)), i) for i in range(n)
                        if rem_of(circle, i) != 0 and circle[i] not in long_arc]
                if not cand:
                    raise ConstructionError(f"deadlock without candidates: {pi} c={c} thr={thr}")
                pos = max(cand)[1]
                long_arc[circle[pos]] = -1 if rem_of(circle, pos) > 0 else 1
                phi = phi_of(circle)
                flips += 1
                continue
            dd, k = step
            if k == 0:
                a, b = head, (head + 1) % n
                circle[a], circle[b] = circle[b], circle[a]
                word.append("X")
                phi = phi_of(circle)
                continue
            if dd != d:
                passes += 1
            d = dd
            head = (head + d) % n
            word.append("L" if d > 0 else "R")
    else:
        raise ValueError(mode)
    # go to c by the shorter arc
    fwd = (c - head) % n
    if fwd <= n - fwd:
        word.append("L" * fwd)
    else:
        word.append("R" * (n - fwd))
    w = "".join(word)
    if apply_word(tuple(pi), w) != identity(n):
        raise ConstructionError(f"word does not sort: {pi} c={c} dir={dirn} thr={thr} {w}")
    nx = w.count("X")
    return {"word": w, "len": len(w), "N_X": nx, "N_rot": len(w) - nx,
            "reduced_len": len(freely_reduce(w)), "passes": passes,
            "c": c, "dir": dirn, "thr": thr, "mode": mode,
            "horizon": horizon if mode == "horizon" else None,
            "weight": weight if mode == "horizon" else None,
            "fallback": fallback if mode == "horizon" else None,
            "long_arc_flips": flips}


def best_sweep(pi, thrs=(1, 2), modes=("pass", "nearest"), horizon=3, weight=1.0, cs=None,
               fallback=False):
    """Shortest sweep word over c, dir, thr, mode (ties: first found).
    cs: iterable of shifts to try (default all)."""
    n = len(pi)
    best = None
    for mode in modes:
        for c in (range(n) if cs is None else cs):
            for dirn in (1, -1):
                for thr in thrs:
                    try:
                        r = sweep_word(pi, c, dirn, thr, mode, horizon=horizon, weight=weight,
                                       fallback=fallback)
                    except ConstructionError:
                        continue
                    if best is None or r["len"] < best["len"]:
                        best = r
    if best is None:
        raise ConstructionError(f"no sweep word for {pi}")
    return best


def _selftest():
    import itertools
    from moves import sigma
    for n in range(4, 8):
        for pi in itertools.permutations(range(n)):
            best_sweep(pi)
    for n in (20, 51, 100):
        best_sweep(sigma(n))
    print("candidate_sweep selftest ok", CANDIDATE_VERSION, CORE_VERSION)


if __name__ == "__main__":
    _selftest()
