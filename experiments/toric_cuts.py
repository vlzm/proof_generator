"""Toric double cuts (H13-I): exact identities and the no-go for averaging families.

Setting (docs/notes/h13_line_model.md).  For a permutation pi of Z_n a double cut
(p, s) produces the line w_j = (pi[(p + j) mod n] - s) mod n, j = 0..n-1, and
I(pi) = min over the n^2 cuts of inv(w).  H13-I claims I(pi) <= floor((n-1)^2/4).

This module is the discovery tool of session 9.  It

  * computes the full n x n table inv(p, s) in O(n^2) per permutation
    (incremental update of the value cut) and I(pi);
  * checks the identities proved in docs/proofs/C37_toric_reflection.md:
      - per-pair cut count:  #{cuts inverting a pair} = n*dx + n*dy - 2*dx*dy;
      - weighted covering formula for a fixed position cut:
          inv(p, s) = sum_i lambda_i * [s in (y_i, y_{i+1}]] - Omega(p),
          lambda_i = (i+1)(n-1-i), Omega(p) independent of s;
      - toric invariance of the triple statistic T(pi) and I(pi) >= T/(n-2);
  * evaluates every candidate *averaging* family of cuts (uniform, anchored at
    the points of pi, all slope families s = a*p + tau) and reports the worst
    permutation for each: all of them fail from small n on, and on affine inputs
    they fail for a structural reason (inv depends only on a*p - s).

Usage: python3 experiments/toric_cuts.py [--nmax 8] [--exhaustive-max 7]
Output: data/runs/toric_cuts/report.md, report.json.
"""

import argparse
import itertools
import json
import os
import random
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

VERSION = "toric_cuts-1.0"
OUT = os.path.join(ROOT, "data", "runs", "toric_cuts")


# ---------------------------------------------------------------- basic tools

def bound(n):
    """floor((n-1)^2/4): the claimed maximum of I(pi)."""
    return ((n - 1) ** 2) // 4


def inv_of(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def line(pi, n, p, s):
    return [(pi[(p + j) % n] - s) % n for j in range(n)]


def inv_table(pi, n):
    """T[p][s] = inv of the line of the cut (p, s).  O(n^2) after O(n^2) start."""
    T = [[0] * n for _ in range(n)]
    for p in range(n):
        vals = [pi[(p + j) % n] for j in range(n)]
        cur = inv_of(vals)
        T[p][0] = cur
        pos = [0] * n
        for j, v in enumerate(vals):
            pos[v] = j
        for s in range(n):
            # shifting the value cut by 1 turns the element of value s (line
            # position m) from the minimum into the maximum
            cur += n - 1 - 2 * pos[s]
            T[p][(s + 1) % n] = cur
    return T


def I_of(pi, n):
    return min(min(row) for row in inv_table(pi, n))


def triples_T(pi, n):
    """Number of triples on which the position and value cyclic orders disagree."""
    T = 0
    for a, b, c in itertools.combinations(range(n), 3):
        seq = (pi[a], pi[b], pi[c])
        k = ((seq[0] > seq[1]) + (seq[0] > seq[2]) + (seq[1] > seq[2]))
        if k % 2 == 1:
            T += 1
    return T


# ------------------------------------------------------------------ identities

def check_pair_formula(n, pi):
    """#cuts (p,s) making the pair inverted = n*dx + n*dy - 2*dx*dy."""
    pts = [(i, pi[i]) for i in range(n)]
    for (x1, y1), (x2, y2) in itertools.combinations(pts, 2):
        dx = (x2 - x1) % n
        dy = (y2 - y1) % n
        cnt = 0
        for p in range(n):
            for s in range(n):
                if (((x1 - p) % n < (x2 - p) % n) !=
                        ((y1 - s) % n < (y2 - s) % n)):
                    cnt += 1
        if cnt != n * dx + n * dy - 2 * dx * dy:
            return False
    return True


def check_weight_formula(n, pi):
    """inv(p,s) = sum_i lambda_i [s in (y_i, y_{i+1}]] - Omega(p)."""
    lam = [(i + 1) * (n - 1 - i) for i in range(n - 1)]
    T = inv_table(pi, n)
    for p in range(n):
        y = [pi[(p + j) % n] for j in range(n)]
        Y = [y[0]]
        for i in range(n - 1):
            Y.append(Y[-1] + (y[i + 1] - y[i]) % n)
        Om = 0
        for j in range(n):
            for k in range(j + 1, n):
                Om += (Y[k] - Y[j] - (y[k] - y[j]) % n) // n
        for s in range(n):
            rhs = -Om
            for i in range(n - 1):
                d = (s - y[i]) % n
                if d != 0 and d <= (y[i + 1] - y[i]) % n:
                    rhs += lam[i]
            if rhs != T[p][s]:
                return False
    return True


def check_triples(n, pi, rng):
    """T is a toric invariant and (n-2) * I(pi) >= T(pi)."""
    T = triples_T(pi, n)
    for _ in range(3):
        a, b = rng.randrange(n), rng.randrange(n)
        pi2 = tuple((pi[(i + a) % n] + b) % n for i in range(n))
        if triples_T(pi2, n) != T:
            return False, "T not toric-invariant"
    if (n - 2) * I_of(pi, n) < T:
        return False, "I < T/(n-2)"
    return True, ""


# ------------------------------------------------------- averaging families

def family_minima(pi, n):
    """Minimal *average* inversion count over several families of n cuts.

    Returns a dict family -> min over the family parameter of the sum over the
    n cuts of the family (compare with n * floor((n-1)^2/4)).
    """
    T = inv_table(pi, n)
    res = {}
    res["uniform"] = sum(sum(row) for row in T) * n // (n * n)  # sum over n^2 cuts / n
    res["diag"] = min(sum(T[p][(p + t) % n] for p in range(n)) for t in range(n))
    res["antidiag"] = min(sum(T[p][(t - p) % n] for p in range(n)) for t in range(n))
    res["slope"] = min(min(sum(T[p][(a * p + t) % n] for p in range(n))
                           for t in range(n)) for a in range(n))
    # anchored: cut (x_A + alpha, y_A + beta) over the n points A of pi
    best = None
    for al in range(n):
        for be in range(n):
            tot = sum(T[(i + al) % n][(pi[i] + be) % n] for i in range(n))
            best = tot if best is None else min(best, tot)
    res["anchored"] = best
    return res


def affine_perms(n):
    out = []
    for a in range(1, n):
        if all((a * i) % n != 0 for i in range(1, n)):
            for b in range(n):
                out.append(tuple((a * i + b) % n for i in range(n)))
    return out


# ------------------------------------------------------------------- driver

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--exhaustive-max", type=int, default=7)
    ap.add_argument("--sample", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=20260915)
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    rng = random.Random(args.seed)
    lines = []
    out = {"version": VERSION, "args": vars(args), "n": {}}

    def log(msg):
        print(msg)
        lines.append(msg)

    log(f"== {VERSION} args={vars(args)}")
    ok = True
    for n in range(4, args.nmax + 1):
        t0 = time.time()
        exhaustive = n <= args.exhaustive_max
        if exhaustive:
            perms = list(itertools.permutations(range(n)))
            cover = f"exhaustive ({len(perms)} perms)"
        else:
            perms = [tuple(rng.sample(range(n), n)) for _ in range(args.sample)]
            perms += affine_perms(n)
            perms += [tuple((h - i) % n for i in range(n)) for h in range(n)]
            cover = f"sampled ({len(perms)} perms: random + affine + reflections)"

        # identities (on a few permutations: the check is O(n^4))
        idpi = [perms[rng.randrange(len(perms))] for _ in range(5)]
        idpi.append(tuple((-i) % n for i in range(n)))
        id_ok = all(check_pair_formula(n, pi) for pi in idpi)
        wf_ok = all(check_weight_formula(n, pi) for pi in idpi)
        tr_ok = all(check_triples(n, pi, rng)[0] for pi in idpi)

        maxI = 0
        argmaxI = None
        worst = {k: (-1, None) for k in
                 ("uniform", "diag", "antidiag", "slope", "anchored")}
        for pi in perms:
            I = I_of(pi, n)
            if I > maxI:
                maxI, argmaxI = I, pi
            fam = family_minima(pi, n)
            for k, v in fam.items():
                if v > worst[k][0]:
                    worst[k] = (v, pi)
        M = bound(n)
        tgt = n * M
        log(f"n={n} {cover}: max I = {maxI} (floor((n-1)^2/4) = {M}) at {argmaxI}; "
            f"identities pair={id_ok} weight={wf_ok} triples={tr_ok}; {time.time()-t0:.1f} s")
        if maxI > M or not (id_ok and wf_ok and tr_ok):
            ok = False
            log(f"  FAIL at n={n}")
        for k in ("uniform", "diag", "antidiag", "slope", "anchored"):
            v, pi = worst[k]
            verdict = "within bound" if v <= tgt else "FAILS the bound"
            log(f"  family {k:9s}: max_pi (sum over the best family) = {v} "
                f"vs n*floor((n-1)^2/4) = {tgt} -> {verdict}; worst pi = {pi}")
        out["n"][n] = {
            "coverage": cover, "max_I": maxI, "argmax_I": list(argmaxI),
            "bound": M, "target_family_sum": tgt,
            "identities": {"pair": id_ok, "weight": wf_ok, "triples": tr_ok},
            "families": {k: {"max_min_sum": worst[k][0], "worst_pi": list(worst[k][1])}
                         for k in worst},
        }

    # affine no-go: inv(p, s) depends only on a*p - s
    log("")
    log("affine no-go (inv depends only on r = a*p + b - s):")
    for n in range(4, args.nmax + 1):
        bad = False
        for pi in affine_perms(n):
            a = (pi[1] - pi[0]) % n
            b = pi[0]
            T = inv_table(pi, n)
            byr = {}
            for p in range(n):
                for s in range(n):
                    r = (a * p + b - s) % n
                    byr.setdefault(r, set()).add(T[p][s])
            if any(len(v) > 1 for v in byr.values()) or \
               any(sum(1 for p in range(n) for s in range(n)
                       if (a * p + b - s) % n == r) != n for r in range(n)):
                bad = True
        log(f"  n={n}: {'FAIL' if bad else 'confirmed'} "
            f"({len(affine_perms(n))} affine inputs; every cut family meeting each "
            f"residue r equally averages to mean_r inv)")
        if bad:
            ok = False
    # the reflection mean
    for n in range(4, args.nmax + 1):
        pi = tuple((-i) % n for i in range(n))
        T = inv_table(pi, n)
        tot = sum(sum(row) for row in T)
        mean = tot / (n * n)
        log(f"  reflection n={n}: mean inv over all cuts = {mean:.3f} "
            f"= C(n,2) - (n^2-1)/6 = {n*(n-1)/2 - (n*n-1)/6:.3f}; "
            f"min = {bound(n)}")

    log(f"verdict: {'PASS' if ok else 'FAIL'}")
    out["verdict"] = "PASS" if ok else "FAIL"
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(out, f, indent=1)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# toric_cuts — двойные разрезы: тождества и запрет усреднений\n\n")
        f.write(f"Команда: `python3 experiments/toric_cuts.py --nmax {args.nmax} "
                f"--exhaustive-max {args.exhaustive_max} --sample {args.sample} "
                f"--seed {args.seed}`. Версия: {VERSION}.\n\n```text\n")
        f.write("\n".join(lines))
        f.write("\n```\n")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
