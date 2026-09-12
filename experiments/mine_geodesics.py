"""Geodesics from id_n to sigma_n (PLAN.md §9.6, claim C13).

Builds the geodesic DAG (states s with d_id(s) + d_sigma(s) = D_n), then a
forward DP counting geodesics by their number of X moves. Verifies C13:
min N_X over geodesics equals floor((n-1)^2/4) with N_rot = floor(n^2/4),
and for 5 <= n the observed splits are reported exactly.
"""

import json
import os
import sys
from collections import deque

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "oracle"))
from moves import apply_move, identity, sigma  # noqa: E402


def bfs_from(start):
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


def mine(n):
    idn, sg = identity(n), sigma(n)
    d_id = bfs_from(idn)
    d_sg = bfs_from(sg)
    D = d_id[sg]

    on_dag = {p for p, d in d_id.items() if d + d_sg[p] == D}
    layers = [[] for _ in range(D + 1)]
    for p in on_dag:
        layers[d_id[p]].append(p)

    # counts[p] = {n_x: number of geodesic prefixes id -> p with n_x swaps}
    counts = {idn: {0: 1}}
    for d in range(D):
        for p in layers[d]:
            cp = counts.get(p)
            if cp is None:
                continue
            for m in ("L", "R", "X"):
                t = apply_move(p, m)
                if t in on_dag and d_id[t] == d + 1:
                    ct = counts.setdefault(t, {})
                    dx = 1 if m == "X" else 0
                    for nx, c in cp.items():
                        ct[nx + dx] = ct.get(nx + dx, 0) + c

    final = counts[sg]
    total = sum(final.values())
    splits = sorted((nx, D - nx, c) for nx, c in final.items())
    min_nx = min(final)
    expect_nx = (n - 1) ** 2 // 4
    expect_rot = n * n // 4
    c13_min = (min_nx == expect_nx and D - min_nx == expect_rot)
    c13_all = (len(final) == 1) if n >= 5 else None
    return {
        "n": n, "D_n": D, "num_geodesics": total,
        "splits_NX_Nrot_count": splits,
        "min_N_X": min_nx,
        "expected_N_X": expect_nx, "expected_N_rot": expect_rot,
        "C13_min_split_ok": c13_min,
        "C13_all_geodesics_same_split_n_ge_5": c13_all,
        "dag_size": len(on_dag),
    }


def main():
    out_dir = os.path.join(os.path.dirname(__file__), "..", "data", "runs",
                           "geodesics")
    os.makedirs(out_dir, exist_ok=True)
    for n in [int(a) for a in sys.argv[1:]]:
        res = mine(n)
        with open(os.path.join(out_dir, f"geodesics_n{n}.json"), "w") as f:
            json.dump(res, f, indent=1)
        print(f"n={n}: D={res['D_n']}, geodesics={res['num_geodesics']}, "
              f"splits (N_X,N_rot,count)={res['splits_NX_Nrot_count']}, "
              f"C13 min split ok: {res['C13_min_split_ok']}, "
              f"all same split (n>=5): "
              f"{res['C13_all_geodesics_same_split_n_ge_5']}", flush=True)


if __name__ == "__main__":
    main()
