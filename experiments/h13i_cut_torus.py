"""h13i_cut_torus-1.0 -- structure of the cut table of the line model (H13-I).

For a permutation pi of Z_n and a double cut (a, b) the line is

    w_j = (pi[(j + a) % n] - b) % n,      j = 0 .. n-1,

and G(a, b) = inv(w) is the cut table; I(pi) = min_{a,b} G(a, b).
H13-I claims I(pi) <= floor((n-1)^2/4) for every pi, n >= 4.

This module checks, exhaustively over all pi for a range of n:

part A (identities, each of them proved in docs/notes/h13_line_model.md 7):
  A1  step recurrences of G against a brute-force recomputation of inv;
  A2  mixed second difference  DD G(a,b) = 2 - 2n*[pi(a) = b];
  A3  closed form  G(a,b) = inv(pi) + A(a) + B(b) + 2ab - 2n*N(a,b);
  A4  discrepancy form  n^2 * G(a,b) = sum(G) - 2*Q(a,b),  Q = n^3 * Dtilde;
  A5  average concordance  sum_{a,b} (B_n - G(a,b)) = sum_{ordered pairs} d*e.

part B (refuted proof routes for H13-I):
  B1  averaging over an anti-diagonal family of cuts {(a, s-a)};
  B2  the relaxation obtained by summing all cut-minimality constraints,
      equivalently sum_i i*pi(i) >= n(n-1)^2/4;
  B3  the one-coordinate strengthening "min_b G(a,b) <= floor((n-1)^2/4)
      for every a";
  B4  the third/second moment bound  max Dtilde >= sum Dt^3 / sum Dt^2.

part C: re-derivation of max_pi I(pi) and of the number of maximizers.

Usage: python3 experiments/h13i_cut_torus.py --nmin 4 --nmax 8 [--out DIR]
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import time
from fractions import Fraction

VERSION = "h13i_cut_torus-1.0"


# ---------------------------------------------------------------- basic tools

def inversions(w) -> int:
    n = len(w)
    c = 0
    for i in range(n):
        wi = w[i]
        for j in range(i + 1, n):
            if wi > w[j]:
                c += 1
    return c


def cut_table(pi):
    """G[a][b] = inv of the line at the double cut (a, b), by the step recurrences.

    G(a+1, b) - G(a, b) = n - 1 - 2 * w_0,   w_0 = (pi[a] - b) mod n
    G(a, b+1) - G(a, b) = n - 1 - 2 * q,     q = (pi^{-1}[b] - a) mod n
    """
    n = len(pi)
    ipi = [0] * n
    for i, v in enumerate(pi):
        ipi[v] = i
    G = [[0] * n for _ in range(n)]
    cur = inversions(list(pi))
    for a in range(n):
        if a > 0:
            cur += (n - 1) - 2 * pi[a - 1]
        G[a][0] = cur
        v = cur
        row = G[a]
        for b in range(1, n):
            v += (n - 1) - 2 * ((ipi[b - 1] - a) % n)
            row[b] = v
    return G


def cut_table_brute(pi):
    n = len(pi)
    return [[inversions([(pi[(j + a) % n] - b) % n for j in range(n)])
             for b in range(n)] for a in range(n)]


def box_counts(pi):
    """N[a][b] = #{i < a : pi[i] < b}."""
    n = len(pi)
    N = [[0] * (n + 1) for _ in range(n + 1)]
    for a in range(n):
        for b in range(n + 1):
            N[a + 1][b] = N[a][b] + (1 if pi[a] < b else 0)
    return N


def q_matrix(pi, N):
    """Q[a][b] = n^3 * Dtilde(a,b), integer; Dtilde = double-centred box discrepancy."""
    n = len(pi)
    P = [[n * N[a][b] - a * b for b in range(n)] for a in range(n)]
    rows = [sum(P[a]) for a in range(n)]
    cols = [sum(P[a][b] for a in range(n)) for b in range(n)]
    tot = sum(rows)
    return [[n * n * P[a][b] - n * rows[a] - n * cols[b] + tot for b in range(n)]
            for a in range(n)]


# ---------------------------------------------------------------- experiment

def run(n, brute_limit=5040):
    """Exhaustive pass over S_n. Returns a dict of results for this n."""
    Bn = n * (n - 1) // 2
    target = (n - 1) ** 2 // 4
    res = {
        "n": n, "B_n": Bn, "floor((n-1)^2/4)": target, "floor(n^2/4)": n * n // 4,
        "A1_brute_checked": 0, "A2_ok": True, "A3_ok": True, "A4_ok": True,
        "A5_ok": True,
    }
    do_brute = len(list(range(n))) >= 0 and _factorial(n) <= brute_limit

    max_I = -1
    argmax_I = None
    n_max_I = 0
    sum_I = 0

    b1_worst = Fraction(-1)      # max over pi of min over s of the diagonal average
    b1_arg = None
    b2_max_inv = -1              # max inv under the summed relaxation
    b2_arg = None
    b3_worst = -1                # max over pi of min_b G(0,b)
    b3_arg = None
    b4_worst = Fraction(-1)      # max over pi of the moment bound
    b4_arg = None

    thr2 = Fraction(n * (n - 1) ** 2, 4)

    for pi in itertools.permutations(range(n)):
        G = cut_table(pi)
        if do_brute:
            if cut_table_brute(pi) != G:
                raise AssertionError(f"A1 failed at n={n}, pi={pi}")
            res["A1_brute_checked"] += 1

        # ---- A2: mixed second difference
        for a in range(n):
            for b in range(n):
                D = (G[(a + 1) % n][(b + 1) % n] - G[(a + 1) % n][b]
                     - G[a][(b + 1) % n] + G[a][b])
                if D != 2 - 2 * n * (1 if pi[a] == b else 0):
                    res["A2_ok"] = False
                    raise AssertionError(f"A2 failed at n={n}, pi={pi}, a={a}, b={b}")

        # ---- A3: closed form
        N = box_counts(pi)
        ipi = [0] * n
        for i, v in enumerate(pi):
            ipi[v] = i
        C = G[0][0]
        Aa = [0] * (n + 1)
        Bb = [0] * (n + 1)
        for a in range(n):
            Aa[a + 1] = Aa[a] + (n - 1) - 2 * pi[a]
        for b in range(n):
            Bb[b + 1] = Bb[b] + (n - 1) - 2 * ipi[b]
        for a in range(n):
            for b in range(n):
                if G[a][b] != C + Aa[a] + Bb[b] + 2 * a * b - 2 * n * N[a][b]:
                    res["A3_ok"] = False
                    raise AssertionError(f"A3 failed at n={n}, pi={pi}, a={a}, b={b}")

        # ---- A4: discrepancy form
        Q = q_matrix(pi, N)
        tot = sum(sum(r) for r in G)
        for a in range(n):
            for b in range(n):
                if n * n * G[a][b] != tot - 2 * Q[a][b]:
                    res["A4_ok"] = False
                    raise AssertionError(f"A4 failed at n={n}, pi={pi}, a={a}, b={b}")

        # ---- A5: average concordance = (1/n^2) sum over ordered pairs of d*e
        s = 0
        for p in range(n):
            for q in range(n):
                if p != q:
                    s += ((q - p) % n) * ((pi[q] - pi[p]) % n)
        if n * n * Bn - tot != s:
            res["A5_ok"] = False
            raise AssertionError(f"A5 failed at n={n}, pi={pi}")

        # ---- C: I(pi)
        I = min(min(r) for r in G)
        sum_I += I
        if I > max_I:
            max_I, argmax_I, n_max_I = I, pi, 1
        elif I == max_I:
            n_max_I += 1

        # ---- B1: anti-diagonal averaging
        best = min(sum(G[a][(s0 - a) % n] for a in range(n)) for s0 in range(n))
        val = Fraction(best, n)
        if val > b1_worst:
            b1_worst, b1_arg = val, pi

        # ---- B2: the summed relaxation
        if Fraction(sum(i * pi[i] for i in range(n))) >= thr2:
            iv = G[0][0]
            if iv > b2_max_inv:
                b2_max_inv, b2_arg = iv, pi

        # ---- B3: one-coordinate strengthening (a = 0 fixed by rotation symmetry)
        v = min(G[0])
        if v > b3_worst:
            b3_worst, b3_arg = v, pi

        # ---- B4: moment bound  max Dtilde >= sum Dt^3 / sum Dt^2
        s2 = 0
        s3 = 0
        for a in range(n):
            for b in range(n):
                q = Q[a][b]
                s2 += q * q
                s3 += q * q * q
        if s2 != 0:
            # G >= avgG - 2n * (s3/s2 scaled back): Q = n^3*Dt, so
            # sum Dt^3 / sum Dt^2 = (s3 / n^9) / (s2 / n^6) = s3 / (n^3 * s2)
            bound = (Fraction(tot, n * n)
                     - 2 * n * Fraction(s3, n ** 3 * s2))
            if bound > b4_worst:
                b4_worst, b4_arg = bound, pi

    res.update({
        "max_I": max_I, "argmax_I": list(argmax_I), "count_max_I": n_max_I,
        "mean_I": float(Fraction(sum_I, _factorial(n))),
        "B1_max_min_diag_avg": str(b1_worst), "B1_value": float(b1_worst),
        "B1_arg": list(b1_arg),
        "B2_max_inv_under_relaxation": b2_max_inv, "B2_arg": list(b2_arg),
        "B3_max_min_over_b": b3_worst, "B3_arg": list(b3_arg),
        "B4_max_moment_bound": str(b4_worst), "B4_value": float(b4_worst),
        "B4_arg": list(b4_arg),
        "B1_refutes": b1_worst > target,
        "B2_refutes": b2_max_inv > target,
        "B3_refutes": b3_worst > target,
        "B4_refutes": b4_worst > target,
    })
    return res


def _factorial(n):
    r = 1
    for i in range(2, n + 1):
        r *= i
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=4)
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--out", default="data/runs/h13i_cut_torus")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    all_res = []
    log = []
    for n in range(args.nmin, args.nmax + 1):
        t0 = time.time()
        r = run(n)
        r["seconds"] = round(time.time() - t0, 1)
        all_res.append(r)
        line = (f"n={n} perms={_factorial(n)} I_max={r['max_I']} "
                f"(target {r['floor((n-1)^2/4)']}, {r['count_max_I']} maximizers) "
                f"B1={r['B1_max_min_diag_avg']} B2={r['B2_max_inv_under_relaxation']} "
                f"B3={r['B3_max_min_over_b']} B4={r['B4_value']:.3f} "
                f"[{r['seconds']}s]")
        print(line, flush=True)
        log.append(line)

    with open(os.path.join(args.out, "report.json"), "w") as f:
        json.dump({"version": VERSION, "results": all_res}, f, indent=1)
    with open(os.path.join(args.out, "report.log"), "w") as f:
        f.write(f"{VERSION}\n" + "\n".join(log) + "\n")


if __name__ == "__main__":
    main()
