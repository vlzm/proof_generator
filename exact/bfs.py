"""BFS distance tables for the LRX Cayley graph.

Two implementations:
- bfs_dict: straightforward dict-based BFS over tuples (reference, n <= 9);
- bfs_lehmer: dense bytearray indexed by Lehmer rank (memory-lean, n >= 8).

They must agree exactly on overlapping n (cross-checked in main).

Table file format: dist_n{n}.bin — one byte per permutation, index = Lehmer
rank, value = distance from id_n, 0xFF = unvisited (must not remain in a
complete table). Metadata in dist_n{n}.json.
"""

import hashlib
import json
import os
import sys
import time
from collections import deque

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "oracle"))
from moves import apply_move, identity, CORE_VERSION  # noqa: E402

UNVISITED = 0xFF


def factorials(n):
    f = [1] * (n + 1)
    for i in range(1, n + 1):
        f[i] = f[i - 1] * i
    return f


def rank_perm(p, fact):
    """Lehmer rank of permutation tuple p (0 .. n!-1)."""
    n = len(p)
    elems = list(p)
    r = 0
    for i in range(n):
        x = elems[i]
        smaller = 0
        for j in range(i + 1, n):
            if elems[j] < x:
                smaller += 1
        r += smaller * fact[n - 1 - i]
    return r


def unrank_perm(r, n, fact):
    """Inverse of rank_perm."""
    avail = list(range(n))
    out = []
    for i in range(n - 1, -1, -1):
        f = fact[i]
        idx, r = divmod(r, f)
        out.append(avail.pop(idx))
    return tuple(out)


def bfs_dict(n):
    """Dict BFS; returns {perm_tuple: distance}."""
    start = identity(n)
    dist = {start: 0}
    q = deque([start])
    while q:
        p = q.popleft()
        d = dist[p] + 1
        for m in ("L", "R", "X"):
            np_ = apply_move(p, m)
            if np_ not in dist:
                dist[np_] = d
                q.append(np_)
    return dist


def bfs_lehmer(n, progress=False):
    """Dense BFS; returns bytearray of length n! indexed by Lehmer rank."""
    fact = factorials(n)
    total = fact[n]
    dist = bytearray([UNVISITED]) * 1
    dist = bytearray(b"\xff" * total)
    start = identity(n)
    r0 = rank_perm(start, fact)
    dist[r0] = 0
    frontier = [r0]
    d = 0
    t0 = time.time()
    while frontier:
        d += 1
        if d >= UNVISITED:
            raise OverflowError("distance exceeds byte encoding")
        nxt = []
        for r in frontier:
            p = unrank_perm(r, n, fact)
            # L
            q = p[1:] + p[:1]
            rq = rank_perm(q, fact)
            if dist[rq] == UNVISITED:
                dist[rq] = d
                nxt.append(rq)
            # R
            q = p[-1:] + p[:-1]
            rq = rank_perm(q, fact)
            if dist[rq] == UNVISITED:
                dist[rq] = d
                nxt.append(rq)
            # X
            q = (p[1], p[0]) + p[2:]
            rq = rank_perm(q, fact)
            if dist[rq] == UNVISITED:
                dist[rq] = d
                nxt.append(rq)
        if progress:
            done = sum(1 for _ in ())  # placeholder, avoid O(n!) scan
            print(f"  layer {d}: frontier {len(nxt)}, "
                  f"elapsed {time.time()-t0:.1f}s", flush=True)
        frontier = nxt
    return dist


def save_table(n, dist_bytes, path_dir, extra_meta=None):
    os.makedirs(path_dir, exist_ok=True)
    bin_path = os.path.join(path_dir, f"dist_n{n}.bin")
    with open(bin_path, "wb") as f:
        f.write(dist_bytes)
    sha = hashlib.sha256(dist_bytes).hexdigest()
    meta = {
        "n": n,
        "generators": "L (left cyclic shift), R (right cyclic shift), "
                      "X (swap first two); PROBLEM.md conventions",
        "indexing": "Lehmer rank (rank_perm in exact/bfs.py)",
        "dtype": "uint8, 0xFF = unvisited",
        "core_version": CORE_VERSION,
        "complete": all(b != UNVISITED for b in dist_bytes),
        "max_distance": max(b for b in dist_bytes if b != UNVISITED),
        "sha256": sha,
    }
    if extra_meta:
        meta.update(extra_meta)
    with open(os.path.join(path_dir, f"dist_n{n}.json"), "w") as f:
        json.dump(meta, f, indent=1)
    return meta


def main():
    tables_dir = os.path.join(os.path.dirname(__file__), "..", "data", "tables")
    ns = [int(a) for a in sys.argv[1:]] or list(range(4, 10))
    for n in ns:
        t0 = time.time()
        fact = factorials(n)
        dist = bfs_lehmer(n)
        dt = time.time() - t0
        # cross-check against dict BFS for n <= 8
        if n <= 8:
            dd = bfs_dict(n)
            assert len(dd) == fact[n]
            for p, d in dd.items():
                assert dist[rank_perm(p, fact)] == d
            check = "cross-checked vs bfs_dict"
        else:
            check = "single implementation (certified separately)"
        meta = save_table(n, dist, tables_dir,
                          {"bfs_seconds": round(dt, 2), "cross_check": check})
        print(f"n={n}: |S_n|={fact[n]}, D_n={meta['max_distance']}, "
              f"complete={meta['complete']}, {dt:.1f}s, {check}", flush=True)


if __name__ == "__main__":
    main()
