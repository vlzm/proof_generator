"""Universal line route ("shrinking cocktail", C34): a fixed head path that sorts
EVERY arrangement of a line of n positions with N_X = inv swaps and N_rot = B_n - 1
rotations (plus the walk from the initial head position to edge 0).

Line model (experiments/line_model.py): positions 0..n-1, edges e = {e, e+1} for
e = 0..n-2, the head sees edge e when head == e; a *reduced* word swaps only
inverted pairs (w[e] > w[e+1]).

Route (edges): pass 1 rightwards 0, 1, ..., n-2; pass 2 leftwards n-3, ..., 0;
pass 3 rightwards 1, ..., n-4; pass 4 leftwards n-5, ..., 1; ... (the k-th pass
has n - k comparators; passes alternate direction and each starts next to the
end of the previous one).  Total comparators 1 + 2 + ... + (n-1) = B_n, path
length B_n - 1.  Correctness: pass 1 carries the maximum to position n-1, pass 2
carries the minimum of positions 0..n-2 to position 0, pass 3 the maximum of
positions 1..n-2 to position n-2, and so on (cocktail shaker sort with
shrinking bounds); every swap is at an inverted pair, so exactly inv(w) swaps
are made.  The bound N_rot <= B_n - 1 is tight for the fully reversed line
(inv = B_n forces N_rot >= N_X - 1 = B_n - 1).

As a circle word (cut edge {n-1, 0} never crossed, head starting at circle
position 0 = line index 0): cocktail_word(w) returns the L/R/X word; it is
verified with the reference moves by checks/check_C34.py.
Version line_cocktail-1.0.
"""

import os
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
from moves import apply_word, identity, CORE_VERSION  # noqa: E402

CONSTRUCTION_VERSION = "line_cocktail-1.0"


def cocktail_route(n):
    """Edge sequence of the shrinking cocktail network for a line of n positions."""
    route = []
    lo, hi = 0, n - 2          # current range of edges
    direction = 1
    while lo <= hi:
        if direction == 1:
            route.extend(range(lo, hi + 1))
            hi -= 1
        else:
            route.extend(range(hi, lo - 1, -1))
            lo += 1
        direction = -direction
    return route


def cocktail_word(w, start=0):
    """Circle word (head starts at circle position `start`, moves along the
    line 0..n-1 only, never crosses the cut edge {n-1, 0}) that sorts the line w
    by the shrinking cocktail route, swapping only inverted pairs.
    Returns (word, N_X, N_rot, edges_swapped, head): after the word the circle is
    sorted and the head stands at line index `head` (the array read from the
    head is the identity rotated by `head`); walking back to 0 costs `head`
    more rotations if the head must end at circle position 0."""
    n = len(w)
    w = list(w)
    head = start
    word = []
    swapped = []
    for e in cocktail_route(n):
        # walk inside the line (no wrap-around): head -> e
        if e > head:
            word.append("L" * (e - head))
        elif e < head:
            word.append("R" * (head - e))
        head = e
        if w[e] > w[e + 1]:
            w[e], w[e + 1] = w[e + 1], w[e]
            word.append("X")
            swapped.append(e)
    word = "".join(word)
    nx = word.count("X")
    return word, nx, len(word) - nx, swapped, head


def _selftest():
    import itertools
    for n in range(2, 8):
        B = n * (n - 1) // 2
        for w in itertools.permutations(range(n)):
            word, nx, nrot, _, head = cocktail_word(w)
            inv = sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])
            assert apply_word(tuple(w), word) == tuple((head + i) % n for i in range(n)), (w, word)
            assert apply_word(tuple(w), word + "R" * head) == identity(n), (w, word)
            assert nx == inv and nrot == B - 1, (w, nx, inv, nrot)
    print("line_cocktail selftest ok", CONSTRUCTION_VERSION, CORE_VERSION)


if __name__ == "__main__":
    _selftest()
