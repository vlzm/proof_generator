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
affine inputs.  Both lose to N1 or to the geodesics in the other regime
(session 5, docs/notes/geodesic_structure.md): "nearest" matches the
geodesics' N_X everywhere but overshoots N_rot on reflections; "pass" is
quadratically wasteful on reflections (rule 9 of AGENTS.md).

Mode "horizon" (sweep-1.2, H12(a) unified head policy): a single rule that
needs no branch on the input family.  After each swap (or on an unproductive
edge), score each direction d by "route cost + gain": the total gain >= thr
reachable within `window` steps (default 3), divided by the distance to the
first such edge; if neither direction has anything within the window, fall
back to the same score over the full remaining range (1..n-1) in each
direction (i.e. the total gain reachable on that side divided by the distance
to its first opportunity there).  Ties keep the current direction.  Session 6
(data/runs/head_policy/): this one rule reproduces N1's short oscillations on
reflections (the local window already finds the next productive edge) and the
geodesics' long sweeps on affine inputs (the fallback favours the side with
more remaining work) -- exhaustively over all pi at 4 <= n <= 9 the length
equals the certified diameter exactly (excess 0, stronger than <= B_n), and
on reflections, sigma_n and its rotations, rev_n and the affine family up to
n = 60 it is <= B_n at every case tried (mostly == B_n).  window = 1 is not
enough (n = 8 has an input with excess +2); window = 2 or 3 closes it in
every case tried so far.

`sweep_word(pi, c, dir, thr)` returns the word and rule-15 statistics and
asserts with the reference moves that the word sorts pi.  `best_sweep(pi)`
minimises the (unreduced) length over c, dir, thr in {1, 2} and mode
(window = 3 for "horizon").  This module never reads distance tables.
"""

import os
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))

from moves import apply_word, freely_reduce, identity, delta, CORE_VERSION  # noqa: E402

CANDIDATE_VERSION = "sweep-1.2"


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


def _horizon_score(circle, target, head, d, n, thr, rem, window):
    """Sum of gain >= thr reachable within `window` steps in direction d (capped at
    n - 1), and the distance to the first such edge (None if none in range)."""
    total = 0
    first = None
    limit = min(window, n - 1)
    for k in range(1, limit + 1):
        h = (head + d * k) % n
        a, b = h, (h + 1) % n
        ra, rb = rem(a), rem(b)
        g = (ra > 0) + (rb < 0) - (ra < 0) - (rb > 0)
        if g >= thr:
            total += g
            if first is None:
                first = k
    return total, first


def sweep_word(pi, c, dirn=1, thr=1, mode="pass", max_passes=None, window=3):
    """mode="pass": alternating full passes (sweep-1.0).  mode="nearest": after each
    swap (or when the current edge is unproductive) the head goes towards the
    nearest productive edge, ties keep the current direction (sweep-1.1).
    mode="horizon" (sweep-1.2, H12(a)): after each swap (or on an unproductive
    edge) score each direction d by the total gain >= thr reachable within
    `window` steps divided by the distance to the first such edge ("route cost
    + gain"); if neither direction has a productive edge within the window,
    fall back to the same score over the full remaining range (1..n-1) in each
    direction, i.e. total reachable gain on that side divided by distance to
    its first opportunity.  Ties keep the current direction.  This one rule
    reproduces N1's short oscillations on reflections (local window suffices)
    and the geodesics' two long sweeps on affine inputs (fallback picks the
    richer side) without branching on the input family."""
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
        while phi > 0:
            guard += 1
            if guard > 8 * n * n:
                raise ConstructionError(f"horizon mode does not terminate: {pi} c={c}")
            a, b = head, (head + 1) % n
            ra, rb = rem(a), rem(b)
            gain = (ra > 0) + (rb < 0) - (ra < 0) - (rb > 0)
            if gain >= thr:
                circle[a], circle[b] = circle[b], circle[a]
                word.append("X")
                phi = sum(delta(i, target[circle[i]], n) for i in range(n))
                continue
            gf, kf = _horizon_score(circle, target, head, 1, n, thr, rem, window)
            gb, kb = _horizon_score(circle, target, head, -1, n, thr, rem, window)
            if kf is None and kb is None:
                gf, kf = _horizon_score(circle, target, head, 1, n, thr, rem, n - 1)
                gb, kb = _horizon_score(circle, target, head, -1, n, thr, rem, n - 1)
                if kf is None and kb is None:
                    raise ConstructionError(f"no productive edge: {pi} c={c} thr={thr}")
            sf = -1 if kf is None else gf / kf
            sb = -1 if kb is None else gb / kb
            if sf > sb or (sf == sb and d > 0):
                nd = 1
            else:
                nd = -1
            if nd != d:
                passes += 1
            d = nd
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
            "c": c, "dir": dirn, "thr": thr, "mode": mode, "window": window}


def best_sweep(pi, thrs=(1, 2), modes=("pass", "nearest", "horizon"), windows=(3,)):
    """Shortest sweep word over c, dir, thr, mode (and window, for mode="horizon";
    ties: first found)."""
    n = len(pi)
    best = None
    for mode in modes:
        for c in range(n):
            for dirn in (1, -1):
                for thr in thrs:
                    for window in (windows if mode == "horizon" else (3,)):
                        try:
                            r = sweep_word(pi, c, dirn, thr, mode, window=window)
                        except ConstructionError:
                            continue
                        if best is None or r["len"] < best["len"]:
                            best = r
                        if mode != "horizon":
                            break
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
