"""h13i_identity_probe-1.0

Session 9 (14.09.2026) exploration of H13-I (PLAN §0/§8, torus double-cut
inversions I(pi) <= floor((n-1)^2/4)).

Two checks, both exhaustive for the stated n:

1. Exact identity for the number of concordant pairs under a double cut
   (q, s) (q = position-rotation origin, s = value-rotation origin,
   equivalent to N1's (q, c) via P = q+1, V = q+1-c):

       TotalConcordant(pi; q, s) = q*(n-q) - J_pi(s) + 2*K_pi(q, s)

   where
     - J_pi(s) = TotalConcordant(pi; 0, s): the single-free-rotation
       (value only, position fixed at identity order) concordant-pair count;
     - K_pi(q, s) = sum over position-pairs i1<i2 of
       g(q; i1, i2) * h(s; pi(i1), pi(i2)), with
       f(z; x, y) = 1{(x - z) mod n < (x - y) mod n} (arc indicator),
       g(q; i1, i2) = f(q; i1, i2), h(s; v1, v2) = f(s; v1, v2).

   This refutes the naively hoped-for decoupled form
   q*(n-q) + s*(n-s) - C(n,2) + 2*K(q,s) (checked and found FALSE in an
   earlier draft of this probe): only the position term q*(n-q) is
   pi-independent (because position order is fixed at the identity
   0..n-1), the value term is NOT s*(n-s) in general but J_pi(s), which
   already depends on pi as strongly as the original problem.

2. Single-free-rotation (value only) is insufficient in general: for some
   pi, max_s J_pi(s) < floor(n^2/4), so the double cut (both q and s free)
   is genuinely necessary, not just a convenience. Confirms/extends the
   Lemma C example of h13_line_model.md (n = 7, (0,5,4,3,2,1,6)).

Usage: python3 experiments/h13i_identity_probe.py
"""
import itertools


def f(z, x, y, n):
    return 1 if (x - z) % n < (x - y) % n else 0


def total_concordant(pi, q, s):
    n = len(pi)
    tot = 0
    for i1 in range(n):
        for i2 in range(i1 + 1, n):
            v1, v2 = pi[i1], pi[i2]
            pj1, pj2 = (i1 - q) % n, (i2 - q) % n
            posfirst = i1 if pj1 < pj2 else i2
            vj1, vj2 = (v1 - s) % n, (v2 - s) % n
            valfirst = i1 if vj1 < vj2 else i2
            if posfirst == valfirst:
                tot += 1
    return tot


def J(pi, s):
    n = len(pi)
    return sum(
        1
        for i1 in range(n)
        for i2 in range(i1 + 1, n)
        if (pi[i1] - s) % n < (pi[i2] - s) % n
    )


def K(pi, q, s):
    n = len(pi)
    tot = 0
    for i1 in range(n):
        for i2 in range(i1 + 1, n):
            v1, v2 = pi[i1], pi[i2]
            tot += f(q, i1, i2, n) * f(s, v1, v2, n)
    return tot


def check_identity(nmax):
    for n in range(4, nmax + 1):
        for pi in itertools.permutations(range(n)):
            pi = list(pi)
            for q in range(n):
                Js = {s: J(pi, s) for s in range(n)}
                for s in range(n):
                    tc = total_concordant(pi, q, s)
                    pred = q * (n - q) - Js[s] + 2 * K(pi, q, s)
                    assert tc == pred, (n, pi, q, s, tc, pred)
        print(f"identity: n={n} OK (all {n}! perms, all q,s)")


def check_single_rotation_insufficiency(nmax):
    for n in range(4, nmax + 1):
        thresh = (n * n) // 4
        worst = None
        for pi in itertools.permutations(range(n)):
            m = max(J(list(pi), s) for s in range(n))
            if worst is None or m < worst[0]:
                worst = (m, pi)
        status = "OK (bound holds)" if worst[0] >= thresh else "VIOLATION (double cut needed)"
        print(f"single rotation: n={n} min_pi max_s J = {worst[0]} vs floor(n^2/4)={thresh}: {status} witness={worst[1]}")


if __name__ == "__main__":
    check_identity(7)
    check_single_rotation_insufficiency(9)
