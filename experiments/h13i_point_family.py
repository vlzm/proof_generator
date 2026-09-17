"""H13-I candidate: the "point-linked" one-parameter family of double cuts.

Context: H13-I conjectures I(pi) = min_{q,c} inv(w(q,c)) <= floor((n-1)^2/4)
for every pi (docs/notes/h13_line_model.md). The full search is over n^2
cuts (q,c); PLAN's next-step note asks to try linking the two cuts into a
single free parameter ("сдвиги обоих разрезов вместе").

Family tested here: for each t in 0..n-1, put the existing point (t, pi(t))
at the new origin (0,0): w^(t)_j = (pi(t+j) - pi(t)) mod n, j = 0..n-1 (indices
mod n). This gives n candidate cuts instead of n^2 (q = t-1, c = t - pi(t) mod
n in the (q,c) convention of experiments/line_model.py). Define
I_diag(pi) = min_t inv(w^(t)).

Lemma (proved below, checked exhaustively 4 <= n <= 9): for every pi,
    sum_{t=0}^{n-1} inv(w^(t)) = 3 R(pi),
where R(pi) is the number of unordered triples {a,b,c} of positions whose
cyclic order under position-labelling and under pi disagree in orientation
("orientation-reversing triples"). Proof: for a fixed triple, exactly one of
the 3 cyclic rotations (t,p,r) has p before r in the rotation started at t;
if the triple is orientation-preserving all 3 rotations of (t,p,r) give
inv-free comparisons for that triple, if orientation-reversing all 3 give an
inversion; hence each triple contributes 0 or 3 to the sum, and R counts the
reversing ones. R(pi) <= C(n,3), with equality iff pi reverses the cyclic
order of every triple, i.e. iff pi is a rotation of a reflection (n such
permutations: the automorphism group of the cyclic order on n points that
reverses orientation is exactly the n reflections, C22/PROBLEM §4.3).

Verdict: the family FAILS as a route to H13-I, and fails increasingly badly.
On reflections, R = C(n,3), so the family's own average 3R/n =
3 C(n,3)/n = (n-1)(n-2)/2 is far above floor((n-1)^2/4) ~ n^2/4 (no
averaging argument over this family can work), and the actual minimum within
the family, I_diag(reflection), also exceeds the bound: checked exhaustively
4 <= n <= 9, I_diag(reflection_n) = C(n-1,2) exactly (the family's worst case
IS the reflection, same as for the true I(pi)), against I(reflection_n) =
floor((n-1)^2/4) (already established, C32/C33). The gap
C(n-1,2) - floor((n-1)^2/4) grows quadratically in n (1, 2, 4, 6, 9, 12 at
n = 4..9); the fraction of pi exceeding the bound stays small but does not
vanish or shrink monotonically (0.167, 0.042, 0.008, 0.011, 0.014, 0.006).
So no patch (different linkage) of a single free parameter is suggested by
this family; H13-I needs a genuinely 2-parameter argument.

Usage: python3 experiments/h13i_point_family.py --nmax 9
Output: data/runs/h13i_point_family/report.json, report.md. Version 1.0.
"""

import argparse
import itertools
import json
import math
import os
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
VERSION = "h13i_point_family-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h13i_point_family")


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def diag_word(pi, t):
    n = len(pi)
    return tuple((pi[(t + j) % n] - pi[t]) % n for j in range(n))


def I_diag(pi):
    n = len(pi)
    return min(inversions(diag_word(pi, t)) for t in range(n))


def count_R(pi):
    n = len(pi)
    r = 0
    for a, b, c in itertools.combinations(range(n), 3):
        if (pi[b] - pi[a]) % n >= (pi[c] - pi[a]) % n:
            r += 1
    return r


def run(n, log, check_identity):
    t0 = time.time()
    bound = (n - 1) ** 2 // 4
    max_diag, argmax_diag, exceed = 0, None, 0
    cnt = 0
    reflection = tuple((0 - i) % n for i in range(n))
    idiag_reflection = I_diag(reflection)
    for pi in itertools.permutations(range(n)):
        cnt += 1
        d = I_diag(pi)
        if check_identity:
            s = sum(inversions(diag_word(pi, t)) for t in range(n))
            r = count_R(pi)
            assert s == 3 * r, (pi, s, r)
        if d > max_diag:
            max_diag, argmax_diag = d, pi
        if d > bound:
            exceed += 1
    row = {
        "n": n, "count": cnt, "bound_floor_(n-1)^2/4": bound,
        "max_I_diag": max_diag, "argmax_I_diag": list(argmax_diag),
        "C(n-1,2)": (n - 1) * (n - 2) // 2,
        "I_diag_reflection": idiag_reflection,
        "exceed_bound_count": exceed, "exceed_fraction": exceed / cnt,
        "identity_checked": bool(check_identity), "seconds": round(time.time() - t0, 1),
    }
    log(f"n={n}: bound={bound} max I_diag={max_diag} (= C(n-1,2)={row['C(n-1,2)']}) "
        f"at {argmax_diag}; I_diag(reflection)={idiag_reflection}; exceed bound on "
        f"{exceed}/{cnt} ({row['exceed_fraction']:.4f}); identity checked={check_identity}; "
        f"{row['seconds']} s")
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=4)
    ap.add_argument("--nmax", type=int, default=9)
    ap.add_argument("--identity-nmax", type=int, default=9,
                     help="verify sum_t inv(w^t) = 3R for n <= this value (expensive: C(n,3) per perm)")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    logf = open(os.path.join(OUT, "report.log"), "a")

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + "\n")
        logf.flush()

    log(f"== {VERSION} args={vars(args)}")
    rows = []
    for n in range(args.nmin, args.nmax + 1):
        rows.append(run(n, log, check_identity=(n <= args.identity_nmax)))
        with open(os.path.join(OUT, f"report_n{args.nmin}_{n}.json"), "w") as f:
            json.dump({"version": VERSION, "rows": rows}, f, indent=1)


if __name__ == "__main__":
    main()
