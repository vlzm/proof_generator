"""check_C37.py — independent checker for C37 (toric double-cut identity) and
for the finite data of C38 (H13-I: I(pi) <= floor((n-1)^2/4)).

Everything here is recomputed from the definitions: lines are relabelled cut by
cut and inversions are counted by the double loop, so the check does not reuse
the difference recurrences of experiments/toric_spearman.py or the C enumerator
experiments/toric_inv.c.  Only integer arithmetic is used (ranks are doubled:
s(u) = 2*(u mod n) - (n-1), so S2 = 4*S).

Parts:
  A. Theorem 3 + Lemma 4: 2 n^2 F(a,b) + n S2(a,b) = n^3 (n-1) - 2 sum_{i!=j} d e
     for every cut, every permutation (exhaustive 3 <= n <= 7; samples 8 <= n <= 11).
  B. Lemma 6: R2 = n^2 (n^2 - 1) - 4 * sum_{a,b} F(a,b), and R2 equals S2 summed
     over the n point cuts.
  C. Corollary 7: 2 n * sum_k F(k, pi(k)) = 6 * sum_{a,b} F(a,b) - n^2 (n^2 - 1).
  D. Corollary 8: 2 n^2 * (n(n-1)/4 - Fbar) = n * E(pi), E counted over triples.
  E. Corollary 9: the shifted point-cut average equals the Delta-sum formula.
  F. H13-I exhaustively over toric class representatives (pi(0) = 0), 4 <= n <= NMAX
     (default 8), and on families (reflections, affine, random) up to n = 60.
  G. Criteria coverage (Corollaries 7 and 9) over all classes, 4 <= n <= 8.

Usage: python3 checks/check_C37.py [--nmax 8] [--seed 20260916]
Version check_C37-1.0.
"""
import argparse
import itertools
import random
import time
from math import gcd


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def line(pi, a, b):
    n = len(pi)
    return [(pi[(a + j) % n] - b) % n for j in range(n)]


def F_grid(pi):
    """F(a,b) for all cuts, straight from the definition."""
    n = len(pi)
    return [[inversions(line(pi, a, b)) for b in range(n)] for a in range(n)]


def S2_grid(pi):
    n = len(pi)
    s = [2 * u - (n - 1) for u in range(n)]
    return [[sum(s[(k - a) % n] * s[(pi[k] - b) % n] for k in range(n)) for b in range(n)]
            for a in range(n)]


def sum_de(pi):
    n = len(pi)
    return sum(((j - i) % n) * ((pi[j] - pi[i]) % n)
               for i in range(n) for j in range(n) if i != j)


def orient(x, y, z, n):
    """+1 if, going up from x, we meet y before z; conventions o(x,x,z)=+1, o(x,y,x)=-1."""
    return 1 if (y - x) % n <= (z - x) % n else -1


def triple_E(pi):
    n = len(pi)
    tot = 0
    for i, j, k in itertools.combinations(range(n), 3):
        tot += orient(i, j, k, n) * orient(pi[i], pi[j], pi[k], n)
    return tot


def check_perm(pi, checks=("A", "B", "C", "D", "E")):
    n = len(pi)
    F = F_grid(pi)
    S2 = S2_grid(pi)
    tot = sum(sum(r) for r in F)
    const = n ** 3 * (n - 1) - 2 * sum_de(pi)
    if "A" in checks:
        for a in range(n):
            for b in range(n):
                assert 2 * n * n * F[a][b] + n * S2[a][b] == const, ("A", pi, a, b)
        # Lemma 4: the constant is 2 n^2 Fbar
        assert const == 2 * tot, ("Lemma 4", pi)
    R2 = sum(S2[k][pi[k]] for k in range(n))
    if "B" in checks:
        s = [2 * u - (n - 1) for u in range(n)]
        R2b = sum(s[(k - l) % n] * s[(pi[k] - pi[l]) % n] for k in range(n) for l in range(n))
        assert R2 == R2b, ("B (two forms of R2)", pi)
        assert R2 == n * n * (n * n - 1) - 4 * tot, ("B (Lemma 6)", pi)
    if "C" in checks:
        pt = sum(F[k][pi[k]] for k in range(n))
        assert 2 * n * pt == 6 * tot - n * n * (n * n - 1), ("C", pi)
    if "D" in checks:
        # Fbar = n(n-1)/4 - E/(2n)  <=>  n^2 (n-1) * n - 4 * ... ; integer form:
        # 4 n^2 Fbar = n^3 (n-1) - 2 n E   with  n^2 Fbar = tot
        assert 4 * tot == n ** 3 * (n - 1) - 2 * n * triple_E(pi), ("D", pi)
    if "E" in checks:
        s = [2 * u - (n - 1) for u in range(n)]
        for (sh, t) in ((0, 0), (1, 0), (0, 1), (n // 2, n // 3)):
            lhs = sum(F[(k + sh) % n][(pi[k] + t) % n] for k in range(n))
            dsum = sum(s[(l - k - sh) % n] * s[(pi[l] - pi[k] - t) % n]
                       for k in range(n) for l in range(n))
            # (1/n) sum_k F = Fbar - (2/n^2) * (1/4) * dsum   ->  integer form:
            assert 4 * n * n * lhs == 4 * n * tot - 2 * n * dsum, ("E", pi, sh, t)
    return F, tot


def I_of(pi):
    return min(min(r) for r in F_grid(pi))


def I_via_identity(pi):
    """I(pi) from Theorem 3 (used only for large n spot checks)."""
    n = len(pi)
    s = [2 * u - (n - 1) for u in range(n)]
    best = None
    # S2(a,b) computed incrementally in O(n^2) overall would need the recurrence;
    # here plainly O(n^3), which is fine for the family sizes used.
    for a in range(n):
        for b in range(n):
            S2 = sum(s[(k - a) % n] * s[(pi[k] - b) % n] for k in range(n))
            if best is None or S2 > best:
                best = S2
    const = n ** 3 * (n - 1) - 2 * sum_de(pi)      # = 2 n^2 Fbar
    # F = (const - n*S2) / (2 n^2)
    num = const - n * best
    assert num % (2 * n * n) == 0, ("non-integer F", pi)
    return num // (2 * n * n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8, help="exhaustive range for part F")
    ap.add_argument("--seed", type=int, default=20260916)
    args = ap.parse_args()
    rnd = random.Random(args.seed)
    t0 = time.time()

    # --- A-E: identities, exhaustive for small n, sampled above -------------
    for n in range(3, 8):
        for pi in itertools.permutations(range(n)):
            check_perm(list(pi))
        print(f"[A-E] identities verified for all {n}! permutations, n = {n}", flush=True)
    for n in range(8, 12):
        for _ in range(30):
            pi = list(range(n))
            rnd.shuffle(pi)
            check_perm(pi)
        for h in (0, 1, n // 2):
            check_perm([(h - i) % n for i in range(n)])
        for m in (2, 3, n - 2):
            if gcd(m, n) == 1:
                check_perm([(m * i + 1) % n for i in range(n)])
        print(f"[A-E] identities verified on 30 random + families, n = {n}", flush=True)

    # --- F: H13-I exhaustively over toric classes ---------------------------
    for n in range(4, args.nmax + 1):
        M = ((n - 1) ** 2) // 4
        worst, arg = -1, None
        for tail in itertools.permutations(range(1, n)):
            pi = [0] + list(tail)
            I = I_of(pi)
            assert I <= M, ("H13-I violated", pi, I)
            if I > worst:
                worst, arg = I, pi
        print(f"[F] n = {n}: max I over toric classes = {worst} = M = {M}, argmax {arg}",
              flush=True)

    # --- F': families for larger n -----------------------------------------
    for n in list(range(12, 25)) + [30, 40, 50, 60]:
        M = ((n - 1) ** 2) // 4
        fam = []
        for h in range(n):
            fam.append(("reflection h=%d" % h, [(h - i) % n for i in range(n)]))
        for m in range(1, n):
            if gcd(m, n) == 1:
                fam.append(("affine m=%d" % m, [(m * i) % n for i in range(n)]))
        for r in range(5):
            p = list(range(n))
            rnd.shuffle(p)
            fam.append(("random %d" % r, p))
        mx, mxname = -1, None
        for name, p in fam:
            I = I_via_identity(p)
            if n <= 16:                     # cross-check the identity route
                assert I == I_of(p), ("identity route", name, n)
            assert I <= M, ("H13-I violated", name, n, p, I)
            if I > mx:
                mx, mxname = I, name
        print(f"[F'] n = {n}: {len(fam)} family members, max I = {mx} <= M = {M} ({mxname})",
              flush=True)

    # --- G: coverage of the proved criteria ---------------------------------
    for n in range(4, 9):
        M = ((n - 1) ** 2) // 4
        s = [2 * u - (n - 1) for u in range(n)]
        cov7 = cov9 = tot_cls = 0
        for tail in itertools.permutations(range(1, n)):
            pi = [0] + list(tail)
            tot_cls += 1
            F = F_grid(pi)
            tot = sum(sum(r) for r in F)
            pt = sum(F[k][pi[k]] for k in range(n))
            if pt <= n * M:
                cov7 += 1
            best = min(sum(F[(k + sh) % n][(pi[k] + t) % n] for k in range(n))
                       for sh in range(n) for t in range(n))
            if best <= n * M:
                cov9 += 1
            assert min(min(r) for r in F) <= M
        print(f"[G] n = {n}: criterion 7 covers {cov7}/{tot_cls} classes, "
              f"criterion 9 covers {cov9}/{tot_cls}", flush=True)

    print(f"PASS  ({time.time() - t0:.1f} s)")


if __name__ == "__main__":
    main()
