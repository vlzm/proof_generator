"""Toric cut grid of a permutation: structure identity, criteria, hard classes (H13-I).

Notation (see docs/proofs/C37_toric_cut_identity.md).  For pi: Z_n -> Z_n and a
cut (a, b) the line is w_j = (pi(a + j) - b) mod n; F(a, b) = inv(w);
I(pi) = min F.  For a point k: X_k = (k - a) mod n, Y_k = (pi(k) - b) mod n.
Doubled centred ranks: s(u) = 2 * (u mod n) - (n - 1).

  S2(a, b) = sum_k s(k - a) * s(pi(k) - b)        (= 4 * Spearman covariance)
  Fbar     = average of F over the n^2 cuts
  R2       = sum_{k, k'} s(k - k') * s(pi(k) - pi(k'))  (= S2 summed over point cuts)

Checked here on every toric class (pi(0) = 0), 4 <= n <= NMAX:
  (A) 2 n^2 F(a, b) + n S2(a, b) = n^3 (n - 1) - 2 * sum_{i != j} d_ij e_ij   (all cuts)
  (B) R2 = n^2 (n^2 - 1) - 4 n^2 Fbar
  (C) criterion I(pi) <= 3 Fbar - (n^2 - 1) / 2   (average over the n point cuts)
  (D) criterion with shifted point cuts: I(pi) <= F(s, t)-average over
      {(k + s, pi(k) + t)}, best over (s, t)
  (E) strengthening: min over bijections beta of the average of F(a, beta(a))
      (assignment problem, Hungarian) -- how often it stays <= floor((n-1)^2/4)
Reported: coverage of the criteria, the hard classes left over, extremes.

Usage: python3 experiments/toric_spearman.py --nmax 8 [--out data/runs/toric_inv]
Version toric_spearman-1.0.
"""
import argparse
import itertools
import json
import os
import time


def inversions(w):
    return sum(1 for i in range(len(w)) for j in range(i + 1, len(w)) if w[i] > w[j])


def grid_direct(pi):
    """F(a, b) for all cuts, relabelling every line (independent of the recurrence)."""
    n = len(pi)
    return [[inversions([(pi[(a + j) % n] - b) % n for j in range(n)]) for b in range(n)]
            for a in range(n)]


def grid_fast(pi):
    """F(a, b) by the two difference recurrences, O(n^2)."""
    n = len(pi)
    pinv = [0] * n
    for i, v in enumerate(pi):
        pinv[v] = i
    row = [0] * n
    row[0] = inversions(list(pi))
    for b in range(1, n):
        row[b] = row[b - 1] + n - 1 - 2 * pinv[b - 1]
    F = [row[:]]
    for a in range(1, n):
        v = pi[a - 1]
        row = [row[b] + n - 1 - 2 * ((v - b) % n) for b in range(n)]
        F.append(row[:])
    return F


def sawtooth(n):
    return [2 * u - (n - 1) for u in range(n)]


def S2_grid(pi):
    n = len(pi)
    s = sawtooth(n)
    return [[sum(s[(k - a) % n] * s[(pi[k] - b) % n] for k in range(n)) for b in range(n)]
            for a in range(n)]


def sum_de(pi):
    """sum over ordered pairs i != j of ((j-i) mod n) * ((pi(j)-pi(i)) mod n)."""
    n = len(pi)
    return sum(((j - i) % n) * ((pi[j] - pi[i]) % n)
               for i in range(n) for j in range(n) if i != j)


def R2_of(pi):
    n = len(pi)
    s = sawtooth(n)
    return sum(s[(k - l) % n] * s[(pi[k] - pi[l]) % n] for k in range(n) for l in range(n))


def hungarian(cost):
    """Minimum cost perfect assignment (O(n^3) Jonker-Volgenant style), integer costs."""
    n = len(cost)
    INF = float("inf")
    u = [0] * (n + 1)
    v = [0] * (n + 1)
    p = [0] * (n + 1)
    way = [0] * (n + 1)
    for i in range(1, n + 1):
        p[0] = i
        j0 = 0
        minv = [INF] * (n + 1)
        used = [False] * (n + 1)
        while True:
            used[j0] = True
            i0 = p[j0]
            delta = INF
            j1 = 0
            for j in range(1, n + 1):
                if not used[j]:
                    cur = cost[i0 - 1][j - 1] - u[i0] - v[j]
                    if cur < minv[j]:
                        minv[j] = cur
                        way[j] = j0
                    if minv[j] < delta:
                        delta = minv[j]
                        j1 = j
            for j in range(n + 1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                else:
                    minv[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while j0:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
    total = 0
    for j in range(1, n + 1):
        total += cost[p[j] - 1][j - 1]
    return total


def analyse(pi, with_assignment=True):
    n = len(pi)
    F = grid_fast(pi)
    I = min(min(r) for r in F)
    tot = sum(sum(r) for r in F)          # = n^2 * Fbar
    M = ((n - 1) ** 2) // 4
    out = {"I": I, "n2Fbar": tot, "M": M}
    # (C) point cuts
    pt = sum(F[k][pi[k]] for k in range(n))            # = n * (average over point cuts)
    out["point_cut_avg_x_n"] = pt
    out["crit_point"] = pt <= n * M
    # (D) shifted point cuts
    best = None
    for s in range(n):
        for t in range(n):
            val = sum(F[(k + s) % n][(pi[k] + t) % n] for k in range(n))
            if best is None or val < best:
                best = val
    out["shift_best_x_n"] = best
    out["crit_shift"] = best <= n * M
    if with_assignment:
        out["assign_x_n"] = hungarian(F)
        out["crit_assign"] = out["assign_x_n"] <= n * M
    return out


def check_identities(pi):
    """(A) and (B) exactly, plus grid_fast == grid_direct."""
    n = len(pi)
    F = grid_fast(pi)
    assert F == grid_direct(pi), ("grid recurrence", pi)
    S2 = S2_grid(pi)
    const = n ** 3 * (n - 1) - 2 * sum_de(pi)
    for a in range(n):
        for b in range(n):
            assert 2 * n * n * F[a][b] + n * S2[a][b] == const, ("identity A", pi, a, b)
    tot = sum(sum(r) for r in F)
    assert R2_of(pi) == n * n * (n * n - 1) - 4 * tot, ("identity B", pi)
    # (B) as a statement about point cuts: R2 = sum over point cuts of S2
    assert R2_of(pi) == sum(S2[k][pi[k]] for k in range(n)), ("identity B'", pi)
    return True


def scan(n, with_assignment=True, identities=True):
    M = ((n - 1) ** 2) // 4
    total = 0
    fail = {"point": 0, "shift": 0, "assign": 0}
    hard_shift = []
    maxI = -1
    max_assign, max_assign_pi = -1, None
    t0 = time.time()
    for tail in itertools.permutations(range(1, n)):
        pi = (0,) + tail
        total += 1
        if identities:
            check_identities(pi)
        r = analyse(list(pi), with_assignment)
        maxI = max(maxI, r["I"])
        assert r["I"] <= M, ("H13-I violated", pi)
        if not r["crit_point"]:
            fail["point"] += 1
        if not r["crit_shift"]:
            fail["shift"] += 1
            if len(hard_shift) < 12:
                hard_shift.append({"pi": list(pi), "I": r["I"],
                                   "shift_avg": r["shift_best_x_n"] / n})
        if with_assignment:
            if not r["crit_assign"]:
                fail["assign"] += 1
            if r["assign_x_n"] > max_assign:
                max_assign, max_assign_pi = r["assign_x_n"], list(pi)
    return {"n": n, "M": M, "classes": total, "max_I": maxI,
            "max_assignment_avg": max_assign / n if with_assignment else None,
            "argmax_assignment": max_assign_pi,
            "fail_point_cut": fail["point"], "fail_shifted_point": fail["shift"],
            "fail_assignment": fail["assign"] if with_assignment else None,
            "hard_for_shifted_point": hard_shift,
            "seconds": round(time.time() - t0, 2)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=4)
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--no-assignment", action="store_true")
    ap.add_argument("--no-identities", action="store_true")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    res = []
    for n in range(args.nmin, args.nmax + 1):
        r = scan(n, not args.no_assignment, not args.no_identities)
        print(json.dumps(r)[:2000], flush=True)
        res.append(r)
    if args.out:
        os.makedirs(args.out, exist_ok=True)
        with open(os.path.join(args.out, f"criteria_n{args.nmin}_{args.nmax}.json"), "w") as f:
            json.dump({"version": "toric_spearman-1.0", "results": res}, f, indent=1)


if __name__ == "__main__":
    main()
