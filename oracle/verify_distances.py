"""Independent certificate for a full distance table (PLAN.md §9.3).

Deliberately self-contained: moves, ranking and unranking are re-implemented
here from the specification in PROBLEM.md, without importing oracle/moves.py
or exact/bfs.py, so that a shared bug in the core cannot certify itself.

Certificate conditions for a table h over all of S_n:
 1. every value is a non-negative integer; the unique zero is at id_n;
 2. for every edge {p, q}: |h(p) - h(q)| <= 1;
 3. for every p != id_n there is a neighbour q with h(q) = h(p) - 1.
Together with completeness this proves h(p) = d(p) for all p.

Neighbour set used: pL, pR, pX. The graph is undirected (L and R are mutually
inverse, X is an involution), so scanning these three neighbours of every p
covers every edge in both directions.
"""

import json
import os
import sys


def _factorials(n):
    f = [1] * (n + 1)
    for i in range(1, n + 1):
        f[i] = f[i - 1] * i
    return f


def _rank(p, fact):
    # Independent Lehmer ranking: count-smaller-to-the-right, radix n-1-i.
    n = len(p)
    r = 0
    for i in range(n):
        c = 0
        for j in range(i + 1, n):
            if p[j] < p[i]:
                c += 1
        r += c * fact[n - 1 - i]
    return r


def _unrank(r, n, fact):
    digits = []
    for i in range(n - 1, -1, -1):
        d, r = divmod(r, fact[i])
        digits.append(d)
    avail = list(range(n))
    return [avail.pop(d) for d in digits]


def certify(n, dist, verbose=True):
    fact = _factorials(n)
    total = fact[n]
    if len(dist) != total:
        return False, f"table length {len(dist)} != {n}! = {total}"
    id_rank = _rank(list(range(n)), fact)

    zeros = 0
    for r in range(total):
        h = dist[r]
        if h == 0xFF:
            return False, f"unvisited entry at rank {r}"
        if h == 0:
            zeros += 1
            if r != id_rank:
                return False, f"zero at rank {r} != id rank {id_rank}"
    if zeros != 1:
        return False, f"{zeros} zeros, expected exactly 1"

    for r in range(total):
        h = dist[r]
        p = _unrank(r, n, fact)
        # neighbours: pL, pR, pX built directly from the move definitions
        nb = (
            p[1:] + p[:1],          # L: left cyclic shift
            p[-1:] + p[:-1],        # R: right cyclic shift
            [p[1], p[0]] + p[2:],   # X: swap first two entries
        )
        has_down = False
        for q in nb:
            hq = dist[_rank(q, fact)]
            if abs(h - hq) > 1:
                return False, f"edge violation at rank {r}: {h} vs {hq}"
            if hq == h - 1:
                has_down = True
        if h > 0 and not has_down:
            return False, f"no descending neighbour at rank {r} (h={h})"
        if verbose and r and r % 500000 == 0:
            print(f"  certified {r}/{total}", flush=True)
    return True, "certificate holds"


def certify_corrupted(n):
    """Self-test: the certifier must reject deliberately damaged tables."""
    from array import array
    tables_dir = os.path.join(os.path.dirname(__file__), "..", "data", "tables")
    with open(os.path.join(tables_dir, f"dist_n{n}.bin"), "rb") as f:
        good = bytearray(f.read())
    fails = 0
    for mutate in (
        lambda t: t.__setitem__(len(t) // 2, 0),          # second zero
        lambda t: t.__setitem__(len(t) // 3, 0xFE),       # huge jump
        lambda t: t.__setitem__(1, t[1] + 2),             # edge violation
    ):
        bad = bytearray(good)
        mutate(bad)
        ok, msg = certify(n, bad, verbose=False)
        if ok:
            return False
        fails += 1
    return fails == 3


def main():
    tables_dir = os.path.join(os.path.dirname(__file__), "..", "data", "tables")
    ns = [int(a) for a in sys.argv[1:]]
    for n in ns:
        with open(os.path.join(tables_dir, f"dist_n{n}.bin"), "rb") as f:
            dist = f.read()
        ok, msg = certify(n, dist)
        print(f"n={n}: certificate {'PASS' if ok else 'FAIL'} — {msg}",
              flush=True)
        if not ok:
            sys.exit(1)
        meta_path = os.path.join(tables_dir, f"dist_n{n}.json")
        with open(meta_path) as f:
            meta = json.load(f)
        meta["certified"] = True
        meta["certifier"] = "oracle/verify_distances.py (independent moves)"
        with open(meta_path, "w") as f:
            json.dump(meta, f, indent=1)
    if ns and certify_corrupted(min(ns)):
        print(f"corruption self-test on n={min(ns)}: certifier rejects all "
              f"3 damaged tables — PASS")
    else:
        print("corruption self-test FAILED", flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
