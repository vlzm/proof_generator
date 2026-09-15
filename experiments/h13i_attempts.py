"""Rejected proof strategies for H13-I (session 9) — reproducible refutations.

H13-I: I(pi) = min over the n^2 double cuts of inv(w) <= floor((n-1)^2/4),
where w_j = pi(a + j) - b (mod n).  Each attempt below is a *general proof
scheme*; the script reports the smallest n at which it fails and a witness.

  1. --deletion   k-element deletion induction.  target(n) - target(n-k) is the
     budget; the scheme needs, for every pi, SOME k-element reduction of the
     toric class (delete k points, contract positions and values) with
     I(pi) <= I(pi'') + budget.  k = 1 was already refuted in session 8; k = 2
     is refuted here (n = 8).
  2. --anchored   averaging over the anchored cut family
     F_{al,be} = {(i + al, pi(i) + be) : i in Z_n} (and F_{al,be,ga} with
     pi(i+ga)).  These are the natural n-cut families that are CONSTANT on both
     extremal classes (reflections and rotations), so the target stays tight on
     them; the scheme needs min over (al, be) of the family average <= target.
  3. --compress   compression towards reflections: a move (swap of values at
     cyclically adjacent positions, or of positions of cyclically adjacent
     values) that raises the winding number omega and does not lower I; the
     scheme would give induction on omega with base omega = n-1 (reflections).

Usage: python3 experiments/h13i_attempts.py [--deletion] [--anchored] [--compress]
       [--nmax 8]   (no flag = all three)
Output: data/runs/h13i_toric/attempts.json (merged with previous content).
Version h13i_attempts-1.0.
"""

import argparse
import itertools
import json
import os

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "data", "runs", "h13i_toric")
VERSION = "h13i_attempts-1.0"


def target(n):
    return ((n - 1) ** 2) // 4


def Ival(pi, n):
    """min over all n^2 cuts, via the shift recurrences (Lemma 1)."""
    if n <= 1:
        return 0
    ipi = [0] * n
    for i, v in enumerate(pi):
        ipi[v] = i
    cur = sum(1 for i in range(n) for j in range(i + 1, n) if pi[i] > pi[j])
    best = row = cur
    for a in range(n):
        if a > 0:
            row = row + n - 1 - 2 * pi[a - 1]
        cur = row
        best = min(best, cur)
        for b in range(n - 1):
            cur = cur + n - 1 - 2 * ((ipi[b] - a) % n)
            best = min(best, cur)
    return best


def inv_table(pi, n):
    ipi = [0] * n
    for i, v in enumerate(pi):
        ipi[v] = i
    t = [[0] * n for _ in range(n)]
    t[0][0] = sum(1 for i in range(n) for j in range(i + 1, n) if pi[i] > pi[j])
    for a in range(n - 1):
        t[a + 1][0] = t[a][0] + n - 1 - 2 * pi[a]
    for a in range(n):
        for b in range(n - 1):
            t[a][b + 1] = t[a][b] + n - 1 - 2 * ((ipi[b] - a) % n)
    return t


def winding(pi, n):
    return sum((pi[(i + 1) % n] - pi[i]) % n for i in range(n)) // n


def reduce_pts(pts):
    xs = sorted(p[0] for p in pts)
    ys = sorted(p[1] for p in pts)
    xi = {x: i for i, x in enumerate(xs)}
    yi = {y: i for i, y in enumerate(ys)}
    out = [0] * len(pts)
    for x, y in pts:
        out[xi[x]] = yi[y]
    return out


def deletion(nmax):
    rows = []
    for n in range(5, nmax + 1):
        for k in (1, 2):
            if n - k < 3:
                continue
            budget = target(n) - target(n - k)
            bad = 0
            first = None
            for pi in itertools.permutations(range(n)):
                I = Ival(list(pi), n)
                pts = [(i, pi[i]) for i in range(n)]
                best = -1
                for dele in itertools.combinations(range(n), k):
                    rest = [pts[i] for i in range(n) if i not in dele]
                    best = max(best, Ival(reduce_pts(rest), n - k))
                    if best >= I - budget:
                        break
                if best < I - budget:
                    bad += 1
                    if first is None:
                        first = {"pi": list(pi), "I": I, "best_reduction_I": best,
                                 "needed": I - budget}
            rows.append({"n": n, "k": k, "budget": budget, "violations": bad,
                         "first": first})
            print(f"deletion n={n} k={k} budget={budget}: violations={bad} {first or ''}")
    return rows


def anchored(nmax, gamma=True):
    rows = []
    for n in range(4, nmax + 1):
        worst = (-1, None)
        worst_plain = -1
        for pi in itertools.permutations(range(n)):
            t = inv_table(list(pi), n)
            best = 10 ** 9
            best_plain = 10 ** 9
            gammas = range(n) if gamma else [0]
            for ga in gammas:
                for al in range(n):
                    for be in range(n):
                        s = sum(t[(i + al) % n][(pi[(i + ga) % n] + be) % n] for i in range(n))
                        best = min(best, s)
                        if ga == 0:
                            best_plain = min(best_plain, s)
            if best > worst[0]:
                worst = (best, list(pi))
            worst_plain = max(worst_plain, best_plain)
        row = {"n": n, "target": target(n),
               "max_min_family_sum": worst[0], "witness": worst[1],
               "max_min_family_sum_gamma0": worst_plain,
               "needed_at_most": n * target(n),
               "holds": worst[0] <= n * target(n)}
        rows.append(row)
        print(f"anchored n={n}: max_pi min_(al,be,ga) avg = {worst[0]}/{n} "
              f"= {worst[0]/n:.3f}, target={target(n)} -> "
              f"{'OK' if row['holds'] else 'FAIL'} (witness {worst[1]})")
    return rows


def compress(nmax):
    rows = []
    for n in range(4, nmax + 1):
        bad = 0
        first = None
        for pi in itertools.permutations(range(n)):
            pi = list(pi)
            w = winding(pi, n)
            if w == n - 1:
                continue
            I = Ival(pi, n)
            ok = False
            cand = []
            for i in range(n):
                q = list(pi)
                q[i], q[(i + 1) % n] = q[(i + 1) % n], q[i]
                cand.append(q)
            ip = [0] * n
            for i, v in enumerate(pi):
                ip[v] = i
            for v in range(n):
                q = list(pi)
                q[ip[v]], q[ip[(v + 1) % n]] = q[ip[(v + 1) % n]], q[ip[v]]
                cand.append(q)
            for q in cand:
                if winding(q, n) > w and Ival(q, n) >= I:
                    ok = True
                    break
            if not ok:
                bad += 1
                if first is None:
                    first = {"pi": pi, "omega": w, "I": I}
        rows.append({"n": n, "violations": bad, "first": first})
        print(f"compress n={n}: pi with no omega-increasing, I-non-decreasing move: {bad} "
              f"{first or ''}")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--deletion", action="store_true")
    ap.add_argument("--anchored", action="store_true")
    ap.add_argument("--compress", action="store_true")
    ap.add_argument("--nmax", type=int, default=8)
    a = ap.parse_args()
    run_all = not (a.deletion or a.anchored or a.compress)
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, "attempts.json")
    data = json.load(open(path)) if os.path.exists(path) else {}
    data["version"] = VERSION
    if a.deletion or run_all:
        data["deletion"] = deletion(a.nmax)
    if a.anchored or run_all:
        data["anchored"] = anchored(min(a.nmax, 7))
    if a.compress or run_all:
        data["compress"] = compress(min(a.nmax, 7))
    with open(path, "w") as f:
        json.dump(data, f, indent=1)
    print("written", path)


if __name__ == "__main__":
    main()
