"""Mine worst-case data from certified distance tables (PLAN.md §9.5).

Checks, per n: D_n and its argmax set; layers D_n-1 and D_n-2; distance
profile; d(rev_n); uniqueness of sigma_n; whether the D_n-1 layer is exactly
the three neighbours of sigma_n; the two distance-preserving symmetries of
PROBLEM.md §4.3 (grouping only, inversion is NOT an automorphism).
"""

import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "oracle"))
from moves import apply_move, identity, sigma, rev  # noqa: E402
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "exact"))
from bfs import factorials, rank_perm, unrank_perm  # noqa: E402


def load(n):
    tables_dir = os.path.join(os.path.dirname(__file__), "..", "data", "tables")
    with open(os.path.join(tables_dir, f"dist_n{n}.bin"), "rb") as f:
        dist = f.read()
    with open(os.path.join(tables_dir, f"dist_n{n}.json")) as f:
        meta = json.load(f)
    assert meta.get("certified"), f"table n={n} not certified"
    return dist, meta


def mine(n):
    dist, _ = load(n)
    fact = factorials(n)
    D = max(dist)
    profile = [0] * (D + 1)
    for b in dist:
        profile[b] += 1
    argmax = [unrank_perm(r, n, fact) for r, b in enumerate(dist) if b == D]
    layer1 = [unrank_perm(r, n, fact) for r, b in enumerate(dist) if b == D - 1]

    s = sigma(n)
    rv = rev(n)
    d_sigma = dist[rank_perm(s, fact)]
    d_rev = dist[rank_perm(rv, fact)]
    sigma_neighbours = sorted(apply_move(s, m) for m in ("L", "R", "X"))

    # PROBLEM §4.3 invariances: d(p) = d(p^{-1}) and d(sigma p sigma)
    inv_ok = conj_ok = True
    for r, b in enumerate(dist):
        p = unrank_perm(r, n, fact)
        pinv = [0] * n
        for i, v in enumerate(p):
            pinv[v] = i
        if dist[rank_perm(tuple(pinv), fact)] != b:
            inv_ok = False
            break
        conj = tuple(s[p[s[i]]] for i in range(n))
        if dist[rank_perm(conj, fact)] != b:
            conj_ok = False
            break

    result = {
        "n": n,
        "D_n": D,
        "B_n": n * (n - 1) // 2,
        "profile": profile,
        "num_at_D": profile[D],
        "num_at_D_minus_1": profile[D - 1],
        "argmax_is_unique_sigma": argmax == [s],
        "layer_Dm1_is_3_neighbours_of_sigma":
            sorted(layer1) == sigma_neighbours,
        "d_sigma": d_sigma,
        "d_rev": d_rev,
        "d_rev_eq_D_minus_2": d_rev == D - 2,
        "inversion_preserves_d": inv_ok,
        "sigma_conjugation_preserves_d": conj_ok,
    }
    return result


def main():
    out_dir = os.path.join(os.path.dirname(__file__), "..", "data", "worst")
    os.makedirs(out_dir, exist_ok=True)
    for n in [int(a) for a in sys.argv[1:]]:
        res = mine(n)
        with open(os.path.join(out_dir, f"worst_n{n}.json"), "w") as f:
            json.dump(res, f, indent=1)
        print(f"n={n}: D_n={res['D_n']} (B_n={res['B_n']}), "
              f"unique sigma argmax: {res['argmax_is_unique_sigma']}, "
              f"layer D-1 = 3 nbrs of sigma: "
              f"{res['layer_Dm1_is_3_neighbours_of_sigma']}, "
              f"d(rev)={res['d_rev']} (=D-2: {res['d_rev_eq_D_minus_2']}), "
              f"inv-sym: {res['inversion_preserves_d']}, "
              f"conj-sym: {res['sigma_conjugation_preserves_d']}", flush=True)


if __name__ == "__main__":
    main()
