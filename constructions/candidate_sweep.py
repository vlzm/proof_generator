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

CANDIDATE_VERSION = "sweep-1.1"


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


def sweep_word(pi, c, dirn=1, thr=1, mode="pass", max_passes=None):
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
            "c": c, "dir": dirn, "thr": thr, "mode": mode}


def best_sweep(pi, thrs=(1, 2), modes=("pass", "nearest")):
    """Shortest sweep word over c, dir, thr, mode (ties: first found)."""
    n = len(pi)
    best = None
    for mode in modes:
        for c in range(n):
            for dirn in (1, -1):
                for thr in thrs:
                    try:
                        r = sweep_word(pi, c, dirn, thr, mode)
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
