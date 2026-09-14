"""Candidate construction "frontier": unified head policy for hypothesis H12
(PLAN.md row H12, docs/notes/geodesic_structure.md section 5, "candidate (a)":
choose the next edge by cost of travel to it plus the gain from the swap,
with a lookahead beyond the single nearest edge).

Session 5's sweep-1.1 "nearest" mode already matches the geodesics' N_X (it
does the swap with the best gain >= thr whenever one is available, and
between swaps walks to the NEAREST edge where a swap is currently possible).
It is exact (<= B_n) on 4 <= n <= 8, but wastes about 0.8*n rotations on
reflections for growing n (data/runs/sweep_eval/report_perms9_fam40.json)
because a pure "go to whichever side has a closer next swap" rule oscillates
locally instead of committing to a direction the way N1's carrier route does.

Two things were tried here (session 6, H12 step 1):

1. **"farthest" mode** -- choose direction by which side's FARTHEST currently
   productive edge is nearer (the "visit the nearer end of a known range
   first" rule that is optimal for the static problem N1's carrier_route
   solves). Empirically this is IDENTICAL to "nearest": on every case tried
   (3000+ random permutations n = 6..14, all of sigma_n/reflections/the H10
   affine counterexamples up to n = 54, see data/runs/candidate_frontier/
   report.md section 1) the produced word is letter-for-letter the same as
   sweep-1.1's "nearest". The set of productive edges on one side of the head
   is essentially always a single contiguous run in these instances, so its
   nearest and farthest end agree on which side is closer; "farthest" adds no
   information over "nearest" here. This matches PLAN.md's own warning that a
   naive lookahead may "plateau like nearest already did" -- kept as
   `frontier_word`/`best_frontier` for the record, not used further.

2. **"weighted" mode** (the one that works): instead of comparing raw
   distances, compare `score = gain - tw * distance` at the NEAREST
   productive edge on each side (tw = 1.0 fixed), and go towards the higher
   score; ties are broken by a fixed `tie_pref` and the very first step by a
   fixed `init_dir`, both searched over like c and thr. This one extra piece
   of information -- a double swap (gain 2) is worth walking one extra step
   past a single swap (gain 1) -- turns out to reproduce N1's short
   oscillations on reflections (the near single swap right next to the head
   stops looking automatically preferable once a double swap one step
   farther is visible) while preserving the long affine sweeps (there, the
   only nearby option is often a single swap while a double swap sits farther
   along the same direction, so `score` keeps pushing the head onward instead
   of turning back). Numerically (data/runs/candidate_frontier/report.md):
   `<= B_n` on every permutation for 4 <= n <= 9 exhaustively (matching the
   certified tables), and strictly `< B_n` (not just `<=`) on reflections,
   sigma_n, rev_n and the H10 affine counterexamples up to n = 60 -- see the
   report for the exact figures and the one caveat (n = 9 needs `tie_pref`
   and `init_dir` as genuine search parameters, not fixed constants; without
   them a handful of symmetric inputs such as the h = 0 reflection lose
   exactly +1).

`weighted_word(pi, c, thr, tw, tie_pref, init_dir)` builds one word and
asserts with the reference moves that it sorts pi. `best_weighted(pi)`
minimises the length over c, thr in {1, 2}, tie_pref in {+1, -1} and
init_dir in {+1, -1} (tw fixed at 1.0). This module never reads distance
tables.

CANDIDATE_VERSION = "frontier-1.1". Core: oracle-1.0.
"""

import os
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))

from moves import apply_word, freely_reduce, identity, delta, CORE_VERSION  # noqa: E402

CANDIDATE_VERSION = "frontier-1.1"


class ConstructionError(Exception):
    pass


def signed_step(a, b, n):
    d = (b - a) % n
    if 2 * d > n:
        d -= n
    return d


def _gain(ra, rb, thr):
    return (ra > 0) + (rb < 0) - (ra < 0) - (rb > 0)


def _init(pi, c):
    n = len(pi)
    c %= n
    circle = list(pi)
    target = [(e + c) % n for e in range(n)]

    def rem(pos):
        return signed_step(pos, target[circle[pos]], n)
    phi = sum(abs(rem(i)) for i in range(n))
    return n, c, circle, target, rem, phi


def _finish(pi, c, n, circle, head, word):
    fwd = (c - head) % n
    if fwd <= n - fwd:
        word.append("L" * fwd)
    else:
        word.append("R" * (n - fwd))
    w = "".join(word)
    if apply_word(tuple(pi), w) != identity(n):
        raise ConstructionError(f"word does not sort: {pi} c={c} {w}")
    nx = w.count("X")
    return {"word": w, "len": len(w), "N_X": nx, "N_rot": len(w) - nx,
            "reduced_len": len(freely_reduce(w)), "c": c}


# ------------------------------------------------------------- mode "farthest"
# Kept for the record (module docstring point 1): empirically identical to
# sweep-1.1's "nearest" on every case tried, so NOT used by best_weighted.

def _farthest_productive(head, d, n, thr, rem):
    best = None
    for k in range(1, n):
        h = (head + d * k) % n
        a, b = h, (h + 1) % n
        if _gain(rem(a), rem(b), thr) >= thr:
            best = k
    return best


def frontier_word(pi, c, thr=1, max_steps=None):
    """"farthest" mode: direction chosen by the nearer FAR end of the
    productive range. See module docstring point 1 -- empirically a no-op
    over "nearest"."""
    n, c, circle, target, rem, phi = _init(pi, c)
    if max_steps is None:
        max_steps = 8 * n * n
    head = 0
    word = []
    d = 1
    steps = 0
    while phi > 0:
        steps += 1
        if steps > max_steps:
            raise ConstructionError(f"frontier does not terminate: {pi} c={c} thr={thr}")
        a, b = head, (head + 1) % n
        if _gain(rem(a), rem(b), thr) >= thr:
            circle[a], circle[b] = circle[b], circle[a]
            word.append("X")
            phi = sum(abs(rem(i)) for i in range(n))
            continue
        kf = _farthest_productive(head, 1, n, thr, rem)
        kb = _farthest_productive(head, -1, n, thr, rem)
        if kf is None and kb is None:
            raise ConstructionError(f"no productive edge: {pi} c={c} thr={thr}")
        if kb is None or (kf is not None and (kf < kb or (kf == kb and d > 0))):
            nd = 1
        else:
            nd = -1
        d = nd
        head = (head + d) % n
        word.append("L" if d > 0 else "R")
    r = _finish(pi, c, n, circle, head, word)
    r.update({"thr": thr, "mode": "farthest"})
    return r


def best_frontier(pi, thrs=(1, 2)):
    n = len(pi)
    best = None
    for c in range(n):
        for thr in thrs:
            try:
                r = frontier_word(pi, c, thr)
            except ConstructionError:
                continue
            if best is None or r["len"] < best["len"]:
                best = r
    if best is None:
        raise ConstructionError(f"no frontier word for {pi}")
    return best


# --------------------------------------------------------------- mode "weighted"
# The one used by best_weighted (module docstring point 2).

def _nearest_score(head, d, n, thr, rem, tw):
    for k in range(1, n):
        h = (head + d * k) % n
        a, b = h, (h + 1) % n
        g = _gain(rem(a), rem(b), thr)
        if g >= thr:
            return g - tw * k
    return None


def weighted_word(pi, c, thr=1, tw=1.0, tie_pref=1, init_dir=1, max_steps=None):
    """Direction chosen by score = gain - tw*distance at the nearest productive
    edge on each side; ties broken by tie_pref, the first step by init_dir
    (both are genuine search parameters of best_weighted, see docstring)."""
    n, c, circle, target, rem, phi = _init(pi, c)
    if max_steps is None:
        max_steps = 8 * n * n
    head = 0
    word = []
    d = init_dir
    reversals = 0
    steps = 0
    while phi > 0:
        steps += 1
        if steps > max_steps:
            raise ConstructionError(f"weighted does not terminate: {pi} c={c} thr={thr}")
        a, b = head, (head + 1) % n
        if _gain(rem(a), rem(b), thr) >= thr:
            circle[a], circle[b] = circle[b], circle[a]
            word.append("X")
            phi = sum(abs(rem(i)) for i in range(n))
            continue
        sf = _nearest_score(head, 1, n, thr, rem, tw)
        sb = _nearest_score(head, -1, n, thr, rem, tw)
        if sf is None and sb is None:
            raise ConstructionError(f"no productive edge: {pi} c={c} thr={thr}")
        if sf is None:
            nd = -1
        elif sb is None:
            nd = 1
        elif sf > sb:
            nd = 1
        elif sb > sf:
            nd = -1
        else:
            nd = tie_pref
        if nd != d:
            reversals += 1
        d = nd
        head = (head + d) % n
        word.append("L" if d > 0 else "R")
    r = _finish(pi, c, n, circle, head, word)
    r.update({"thr": thr, "tw": tw, "tie_pref": tie_pref, "init_dir": init_dir,
              "reversals": reversals, "mode": "weighted"})
    return r


def best_weighted(pi, thrs=(1, 2), ties=(1, -1), dirs=(1, -1), tw=1.0):
    """Shortest weighted word over c, thr, tie_pref, init_dir (tw fixed)."""
    n = len(pi)
    best = None
    for c in range(n):
        for thr in thrs:
            for tie in ties:
                for idir in dirs:
                    try:
                        r = weighted_word(pi, c, thr, tw, tie, idir)
                    except ConstructionError:
                        continue
                    if best is None or r["len"] < best["len"]:
                        best = r
    if best is None:
        raise ConstructionError(f"no weighted word for {pi}")
    return best


def _selftest():
    import itertools
    from moves import sigma
    for n in range(4, 8):
        for pi in itertools.permutations(range(n)):
            best_weighted(pi)
    for n in (20, 51, 100):
        best_weighted(sigma(n))
    print("candidate_frontier selftest ok", CANDIDATE_VERSION, CORE_VERSION)


if __name__ == "__main__":
    _selftest()
