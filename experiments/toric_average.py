"""toric_average.py — averaging over the n^2 double cuts: exact mean, the
sufficient condition of C37 and the failure rates of the restricted averaging
families (H13-I, session 9).

For every toric class of S_n (one representative per class shift, pi[0] = 0) the
script computes

  * I(pi) = min over the n^2 cuts (a, s) of inv(a, s)  (line w_j = pi(a+j) - s);
  * the exact sum of inv over all cuts, checked against the identity of C37
    Lemma 2:  12 * sum = 3 n^3 (n-1) + n^2 (n-1)(n-2) - 12 Delta(pi),
    Delta(pi) = sum over pairs of ((i-j mod n) + (pi i - pi j mod n) - n)^2;
  * whether the plain average over all n^2 cuts is <= floor((n-1)^2/4)
    (C37 Corollary 4), and the same for the restricted families:
      - lines s = c + lam * a of the cut torus (n^2 + n families of n cuts),
      - anti-diagonals a + s = c (the family that is exactly optimal on
        reflections),
      - "first line element equals t" (n families of n cuts).

Usage: python3 experiments/toric_average.py [--nmax 8] [--out DIR]
Output: DIR/toric_average.json (default data/runs/toric_inversions/).
Version toric_average-1.0.
"""

import argparse
import itertools
import json
import os
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
VERSION = "toric_average-1.0"


def inv_table(pi):
    """inv[a][s] for all n^2 cuts, by the incremental rule (C37 Lemma 7)."""
    n = len(pi)
    pinv = [0] * n
    for i, v in enumerate(pi):
        pinv[v] = i
    base = sum(1 for i in range(n) for j in range(i + 1, n) if pi[i] > pi[j])
    table = []
    cur_a = base
    for a in range(n):
        if a:
            cur_a += n - 1 - 2 * pi[a - 1]
        row = [cur_a]
        cur = cur_a
        for s in range(1, n):
            cur += n - 1 - 2 * ((pinv[s - 1] - a) % n)
            row.append(cur)
        table.append(row)
    return table


def defect(pi):
    n = len(pi)
    return sum(((i - j) % n + (pi[i] - pi[j]) % n - n) ** 2
               for i in range(n) for j in range(i + 1, n))


def scan(n):
    target = (n - 1) ** 2 // 4
    res = {
        "n": n, "target": target, "reps": 0, "max_I": -1, "argmax_I": None,
        "mean_fail": 0, "line_fail": 0, "anti_fail": 0,
        "t_fail": [0] * n, "min_t_sets": [], "seconds": 0.0,
    }
    t0 = time.time()
    bad_t = set()
    for pi in itertools.permutations(range(n)):
        if pi[0] != 0:
            continue
        res["reps"] += 1
        tab = inv_table(pi)
        I = min(min(r) for r in tab)
        if I > res["max_I"]:
            res["max_I"], res["argmax_I"] = I, list(pi)
        total = sum(sum(r) for r in tab)
        assert 12 * total == 3 * n ** 3 * (n - 1) + n * n * (n - 1) * (n - 2) - 12 * defect(pi), \
            f"C37 Lemma 2 failed at n={n}, pi={pi}"
        if total > target * n * n:
            res["mean_fail"] += 1
        best_line = min(min(sum(tab[a][(c + lam * a) % n] for a in range(n))
                            for c in range(n)) for lam in range(n))
        best_line = min(best_line, min(sum(tab[a]) for a in range(n)))
        if best_line > target * n:
            res["line_fail"] += 1
        anti = min(sum(tab[a][(c - a) % n] for a in range(n)) for c in range(n))
        if anti > target * n:
            res["anti_fail"] += 1
        fails = []
        for t in range(n):
            if min(tab[a][(pi[a] - t) % n] for a in range(n)) > target:
                res["t_fail"][t] += 1
                fails.append(t)
        bad_t.add(frozenset(fails))
    for size in (1, 2, 3):
        good = [list(S) for S in itertools.combinations(range(n), size)
                if all(not set(S) <= f for f in bad_t)]
        if good:
            res["min_t_sets"] = {"size": size, "sets": good}
            break
    res["seconds"] = round(time.time() - t0, 1)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--nmin", type=int, default=4)
    ap.add_argument("--out", default=os.path.join(ROOT, "data", "runs", "toric_inversions"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rows = []
    for n in range(args.nmin, args.nmax + 1):
        r = scan(n)
        rows.append(r)
        print(json.dumps(r))
    path = os.path.join(args.out, "toric_average.json")
    old = []
    if os.path.exists(path):
        with open(path) as f:
            old = json.load(f).get("rows", [])
    keep = [r for r in old if not any(r["n"] == x["n"] for x in rows)]
    with open(path, "w") as f:
        json.dump({"version": VERSION, "rows": sorted(keep + rows, key=lambda r: r["n"])}, f, indent=1)
    print("written", path)


if __name__ == "__main__":
    main()
