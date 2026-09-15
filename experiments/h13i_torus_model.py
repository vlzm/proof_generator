"""H13-I torus reformulation and structural identities (session 9).

H13-I conjectures I(pi) = min_{q,c} inv(w_{q,c}) <= floor((n-1)^2/4) for every
permutation pi of Z_n, n >= 4 (see docs/notes/h13_line_model.md, PLAN.md H13
row).  This script independently verifies, by brute force against the same
w_{q,c} definition used in experiments/line_profile.c and
experiments/line_model.py, four structural facts found and used in session 9
while searching for a proof of H13-I.  None of them by itself proves H13-I;
they are recorded here because a future session may be able to complete the
argument, and because rule 7 requires reproducible experiments to be logged.

Coordinates.  Write a = q+1, b = q+1-c (both in Z_n).  Then
    w_j(q,c) = pi((a+j) mod n) - b  (mod n),   j = 0, ..., n-1,
and F(a,b) := inv(w_{q,c}) with q = a-1, c = a-b (all mod n).  This matches
line_profile.c/line_model.py exactly (checked in fact 1 below) and is the
"double cut" (a,b) on the torus Z_n x Z_n described in the session 9 task.

Facts checked (brute force, all pi, small n; see main() for ranges):

 1. F(a,b) as defined above via (a,b) equals inv(w_{q,c}) as defined via
    (q,c) in line_profile.c, for every pi, q, c.  (Coordinate change check.)

 2. Quadrant identity.  For an unordered pair of points {(i1,v1),(i2,v2)}
    of the graph of pi, with position gap l = (i2-i1) mod n and value gap
    m = (v2-v1) mod n (either labelling of the pair; the statement below is
    symmetric under l -> n-l, m -> n-m together), the pair is "discordant"
    (contributes to inv(w_{a,b})) for exactly l*(n-m) + (n-l)*m of the n^2
    cuts (a,b), and "concordant" for the complementary l*m + (n-l)*(n-m).
    Equivalently concordant(pair) - discordant(pair) = (n-2l)*(n-2m), summed
    over cuts.  Checked by direct enumeration of cuts per pair.

 3. Row/column step identities.  Writing F(a,b) for the double-cut inversion
    count,
        F(a+1,b) - F(a,b) = (n-1) - 2*((pi(a) - b) mod n),
        F(a,b+1) - F(a,b) = (n-1) - 2*((pi^{-1}(b) - a) mod n),
    for all a, b in Z (using pi(a mod n) etc.); i.e. moving the position cut
    or the value cut by one step changes inv by an explicit, computable
    amount.  Proved by the standard "move an extremal element from one end
    of a linear arrangement to the other" inversion-count argument (see
    docs/notes/h13_line_model.md, session 9 section, for the derivation).

 4. Closed form (telescoping the two step identities over a rectangle):
    for integers 0 <= a, b <= n (no wraparound),
        F(a,b) = F(a,0) + F(0,b) - F(0,0) + 2*a*b - 2*n*R(a,b),
    where R(a,b) = #{i : 0 <= i < a, 0 <= pi(i) < b} counts points of pi's
    graph in the rectangle [0,a) x [0,b).  Equivalently the discrete mixed
    second difference of F is +2 everywhere except at the n cells whose
    lower-left corner (a,b) satisfies b = pi(a) mod n, where it is
    -2*(n-1).  This is the discrete "Hessian" of F: away from pi's own
    graph, F behaves exactly like the identity permutation's F, which is
    the classical rotation-inversion tent phi(c) = c*(n-c), c = (a-b) mod n.

 5. Arithmetic identity C(n,2) - floor(n^2/4) = floor((n-1)^2/4) (both
    parities of n), i.e. H13-I ("some cut has <= floor((n-1)^2/4)
    inversions") is equivalent to "some cut has >= floor(n^2/4) concordant
    pairs".

None of facts 1-5 by itself proves H13-I; fact 4 does show that F(a,b) is,
away from pi's own n points, identical in shape to the extremal case
(reflections), which is consistent with reflections being the observed
maximizers of I(pi) (C33/C35) but does not by itself bound min_{a,b} F(a,b).
See docs/notes/h13_line_model.md, session 9 section, for what was tried and
why it did not close the gap.

Usage: python3 experiments/h13i_torus_model.py --nmax 8
Version h13i_torus_model-1.0.
"""

import argparse
import itertools
import json
import os
import random
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
VERSION = "h13i_torus_model-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h13i_torus_model")


def inv(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def w_qc(pi, q, c):
    n = len(pi)
    shift = (q + 1 - c) % n
    return [(pi[(q + 1 + j) % n] - shift) % n for j in range(n)]


def F(pi, a, b):
    """Double-cut inversion count in (a, b) = (q+1, q+1-c) coordinates."""
    n = len(pi)
    return inv([(pi[(a + j) % n] - b) % n for j in range(n)])


def check_coordinate_change(n, log):
    for pi in itertools.permutations(range(n)):
        for q in range(n):
            for c in range(n):
                a = (q + 1) % n
                b = (q + 1 - c) % n
                if F(pi, a, b) != inv(w_qc(pi, q, c)):
                    return False, (pi, q, c)
    log(f"n={n}: fact 1 (coordinate change (q,c) <-> (a,b)) OK, all pi, all q, c")
    return True, None


def gap_class(n, i1, v1, i2, v2):
    l = (i2 - i1) % n
    m = (v2 - v1) % n
    return l, m


def check_quadrant_identity(n, log):
    for pi in itertools.permutations(range(n)):
        pts = [(i, pi[i]) for i in range(n)]
        for (i1, v1), (i2, v2) in itertools.combinations(pts, 2):
            l, m = gap_class(n, i1, v1, i2, v2)
            disc = conc = 0
            for a in range(n):
                for b in range(n):
                    j1 = (i1 - a) % n
                    j2 = (i2 - a) % n
                    w1 = (v1 - b) % n
                    w2 = (v2 - b) % n
                    pos_1_first = j1 < j2
                    val_1_first = w1 < w2
                    if pos_1_first == val_1_first:
                        conc += 1
                    else:
                        disc += 1
            expected_disc = l * (n - m) + (n - l) * m
            expected_conc = l * m + (n - l) * (n - m)
            if disc != expected_disc or conc != expected_conc:
                return False, (pi, i1, v1, i2, v2, disc, expected_disc)
    log(f"n={n}: fact 2 (quadrant concordant/discordant counts per pair) OK, all pi, all pairs")
    return True, None


def check_row_col_steps(n, log, trials=200, seed=0):
    rng = random.Random(seed)
    for _ in range(trials):
        pi = list(range(n))
        rng.shuffle(pi)
        pinv = [0] * n
        for i, v in enumerate(pi):
            pinv[v] = i
        a = rng.randrange(0, 3 * n)
        b = rng.randrange(0, 3 * n)
        lhs_a = F(pi, a + 1, b) - F(pi, a, b)
        rhs_a = (n - 1) - 2 * ((pi[a % n] - b) % n)
        lhs_b = F(pi, a, b + 1) - F(pi, a, b)
        rhs_b = (n - 1) - 2 * ((pinv[b % n] - a) % n)
        if lhs_a != rhs_a or lhs_b != rhs_b:
            return False, (pi, a, b, lhs_a, rhs_a, lhs_b, rhs_b)
    log(f"n={n}: fact 3 (row/column step identities) OK, {trials} random (pi, a, b), a,b up to 3n (wraparound exercised)")
    return True, None


def rectangle_count(pi, a, b):
    return sum(1 for i in range(a) if pi[i] < b)


def check_closed_form(n, log, trials=300, seed=1):
    rng = random.Random(seed)
    for _ in range(trials):
        pi = list(range(n))
        rng.shuffle(pi)
        a = rng.randint(0, n)
        b = rng.randint(0, n)
        lhs = F(pi, a, b)
        rhs = F(pi, a, 0) + F(pi, 0, b) - F(pi, 0, 0) + 2 * a * b - 2 * n * rectangle_count(pi, a, b)
        if lhs != rhs:
            return False, (pi, a, b, lhs, rhs)
    log(f"n={n}: fact 4 (closed form F(a,b) = F(a,0)+F(0,b)-F(0,0)+2ab-2nR(a,b)) OK, {trials} random (pi, a, b), 0<=a,b<=n")
    return True, None


def check_arithmetic_identity(nmax, log):
    for n in range(2, nmax + 1):
        c2 = n * (n - 1) // 2
        lhs = c2 - n * n // 4
        rhs = (n - 1) * (n - 1) // 4
        if lhs != rhs:
            return False, n
    log(f"n=2..{nmax}: fact 5 (C(n,2) - floor(n^2/4) = floor((n-1)^2/4)) OK")
    return True, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=4)
    ap.add_argument("--nmax", type=int, default=7)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    logf = open(os.path.join(OUT, "report.log"), "a")

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + "\n")

    t0 = time.time()
    results = {"version": VERSION, "checks": []}
    ok_all = True
    ok, bad = check_arithmetic_identity(30, log)
    ok_all &= ok
    for n in range(args.nmin, args.nmax + 1):
        for name, fn in (("coordinate_change", check_coordinate_change),
                          ("quadrant_identity", check_quadrant_identity)):
            ok, bad = fn(n, log)
            ok_all &= ok
            results["checks"].append({"n": n, "fact": name, "ok": ok, "counterexample": bad})
            if not ok:
                log(f"n={n}: FACT {name} FAILED: {bad}")
        for name, fn in (("row_col_steps", check_row_col_steps),
                          ("closed_form", check_closed_form)):
            ok, bad = fn(n, log)
            ok_all &= ok
            results["checks"].append({"n": n, "fact": name, "ok": ok, "counterexample": bad})
            if not ok:
                log(f"n={n}: FACT {name} FAILED: {bad}")
    results["ok_all"] = ok_all
    results["seconds"] = round(time.time() - t0)
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(results, f, indent=1)
    log(f"all facts 1-5 {'PASS' if ok_all else 'FAIL'}; {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
