"""check_C37.py — independent check of docs/proofs/C37_toric_average.md.

Claim C37 (PROVED part):
  L1  inv(u+1,v) - inv(u,v) = -2*phi(pi(u)-v),  inv(u,v+1) - inv(u,v) = -2*phi(sigma(v)-u)
  L2  inv(u,v) + (2/n)*T(u,v) = Ibar(pi)          (identity on the whole cut torus)
  L3  k_a integer in [1,n-1];  k_a + k_{n-a} = n;  k_a+k_b+k_c <= 2n when a+b+c = 0 mod n
  L4  Ibar = n(n-1)/2 - (1/n) sum_a a*k_a
  L5  sum_a a*k_a >= n(n^2-1)/6
  TH  I(pi) <= Ibar(pi) <= (n-1)(2n-1)/6

Also reported (data for C38, not proved):
  H13-I  I(pi) <= floor((n-1)^2/4)                      exhaustive here for 4<=n<=8
  RA     sum_u min_v inv(u,v) <= n*floor((n-1)^2/4)     exhaustive here for 4<=n<=8

Everything is exact integer arithmetic: phi is used only through
psi(x) = 2*((x mod n)) - (n-1) = 2*phi(x), and the identity is checked in the
integral form  n*(2n*inv(u,v) + T2(u,v)) = 2 * sum_{u,v} inv(u,v),
where T2 = sum_i psi(i-u)*psi(pi(i)-v) = 4*T.

Usage: python3 checks/check_C37.py [--nmax 8] [--sample_nmax 40] [--seed 1]
Version check_C37-1.0.  No dependency on oracle/moves.py: no words are executed
here, only inversion counts of relabelled lines.
"""
import argparse
import itertools
import random
import sys
import time


# ---------------------------------------------------------------- basics

def inv_line(w):
    n = len(w)
    c = 0
    for i in range(n):
        wi = w[i]
        for j in range(i + 1, n):
            if wi > w[j]:
                c += 1
    return c


def line(pi, u, v):
    n = len(pi)
    return [(pi[(u + j) % n] - v) % n for j in range(n)]


def inv_table_direct(pi):
    """inv(u,v) for every cut, computed from scratch (definition)."""
    n = len(pi)
    return [[inv_line(line(pi, u, v)) for v in range(n)] for u in range(n)]


def inv_table_fast(pi, sigma):
    """same table by the step formulas of Lemma 1 (used for the big ladders)."""
    n = len(pi)
    tab = [[0] * n for _ in range(n)]
    tab[0][0] = inv_line(line(pi, 0, 0))
    for v in range(n - 1):
        tab[0][v + 1] = tab[0][v] + n - 1 - 2 * sigma[v]
    for u in range(n - 1):
        for v in range(n):
            tab[u + 1][v] = tab[u][v] + n - 1 - 2 * ((pi[u] - v) % n)
    return tab


def psi(x, n):
    return 2 * (x % n) - (n - 1)


def T2(pi, u, v):
    n = len(pi)
    return sum(psi(i - u, n) * psi(pi[i] - v, n) for i in range(n))


def kvec(pi):
    n = len(pi)
    k = [0] * n
    for a in range(1, n):
        s = sum((pi[(i + a) % n] - pi[i]) % n for i in range(n))
        assert s % n == 0, ("k_a not integral", pi, a)
        k[a] = s // n
    return k


def inverse(pi):
    n = len(pi)
    s = [0] * n
    for i, p in enumerate(pi):
        s[p] = i
    return s


# ---------------------------------------------------------------- checks

def check_perm(pi, full=True):
    """all lemmas for one permutation; returns (I, RA, Ibar_num) with
    Ibar_num = sum over all cuts of inv."""
    n = len(pi)
    sigma = inverse(pi)
    tab = inv_table_fast(pi, sigma)
    if full:
        ref = inv_table_direct(pi)
        assert tab == ref, ("L1 step formulas disagree with the definition", pi)
    tot = sum(sum(r) for r in tab)
    # L2: identity on the whole torus, integral form
    if full:
        for u in range(n):
            for v in range(n):
                assert n * (2 * n * tab[u][v] + T2(pi, u, v)) == 2 * tot, \
                    ("L2 identity fails", pi, u, v)
    else:
        for (u, v) in ((0, 0), (1, 0), (0, 1), (n // 2, n // 3), (n - 1, n - 1)):
            assert n * (2 * n * tab[u % n][v % n] + T2(pi, u % n, v % n)) == 2 * tot, \
                ("L2 identity fails (sampled)", pi, u, v)
    # L3
    k = kvec(pi)
    for a in range(1, n):
        assert 1 <= k[a] <= n - 1, ("L3.1", pi, a, k[a])
        assert k[a] + k[n - a] == n, ("L3.2", pi, a)
    for a in range(1, n):
        for b in range(1, n):
            c = (-a - b) % n
            if c == 0:
                continue
            assert k[a] + k[b] + k[c] <= 2 * n, ("L3.3", pi, a, b, c)
    # L4: Ibar = n(n-1)/2 - (1/n) sum a k_a, checked as n^2*Ibar = tot
    s_ak = sum(a * k[a] for a in range(1, n))
    assert tot * 2 == n * n * n * (n - 1) - 2 * n * s_ak, ("L4", pi, tot, s_ak)
    # L5
    assert 6 * s_ak >= n * (n * n - 1), ("L5", pi, s_ak)
    # theorem: I <= Ibar <= (n-1)(2n-1)/6   (in integers: 6*n^2*Ibar <= n^2*(n-1)(2n-1))
    I = min(min(r) for r in tab)
    assert 6 * tot <= n * n * (n - 1) * (2 * n - 1), ("TH upper", pi, tot)
    assert I * n * n <= tot, ("TH min<=mean", pi)
    RA = sum(min(r) for r in tab)
    return I, RA, tot


def check_cut_convention(n):
    """the (q, c) cuts of docs/notes/h13_line_model.md give the same multiset
    of inv values as the (u, v) cuts used here."""
    for pi in itertools.permutations(range(n)):
        a = sorted(inv_line([(pi[(q + 1 + j) % n] - (q + 1 - c)) % n for j in range(n)])
                   for q in range(n) for c in range(n))
        b = sorted(inv_line(line(pi, u, v)) for u in range(n) for v in range(n))
        assert a == b, ("cut conventions differ", pi)


def check_corollaries(n):
    """Corollaries 3 and 4 of docs/proofs/C37_toric_average.md:
       T = n^2 * (doubly centred corner discrepancy D),
       mixed second difference of inv = 2 - 2n*[pi(u)=v].
    Fractions are avoided: D(u,v) = N(u,v) - u*v/n is used as n*D."""
    for pi in itertools.permutations(range(n)):
        nD = [[n * sum(1 for i in range(u) if pi[i] < v) - u * v for v in range(n)]
              for u in range(n)]
        rows = [sum(nD[u]) for u in range(n)]
        cols = [sum(nD[u][v] for u in range(n)) for v in range(n)]
        tot = sum(rows)
        for u in range(n):
            for v in range(n):
                # 4*T = T2, and n^2*(D - Drow - Dcol + Dbar) = n*(nD) - rows - cols + tot/n
                lhs = n * n * T2(pi, u, v)
                rhs = 4 * n * (n * n * nD[u][v] - n * rows[u] - n * cols[v] + tot)
                assert lhs == rhs, ("corollary 3 (T = n^2 * centred discrepancy)", pi, u, v)
        tab = inv_table_fast(pi, inverse(pi))
        for u in range(n):
            for v in range(n):
                m = (tab[(u + 1) % n][(v + 1) % n] + tab[u][v]
                     - tab[(u + 1) % n][v] - tab[u][(v + 1) % n])
                assert m == 2 - 2 * n * (1 if pi[u] == v else 0), \
                    ("corollary 4 (mixed difference)", pi, u, v, m)


def reflections(n):
    return [tuple((h - i) % n for i in range(n)) for h in range(n)]


def samples(n, seed):
    rnd = random.Random(seed + n)
    out = []
    out.extend(reflections(n)[:3])
    out.append(tuple(range(n)))
    for a in range(1, n):
        if all(a * i % n != 0 for i in range(1, n)):  # gcd(a,n)=1
            out.append(tuple((a * i + 1) % n for i in range(n)))
    for _ in range(12):
        p = list(range(n))
        rnd.shuffle(p)
        out.append(tuple(p))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8, help="exhaustive ladder up to this n")
    ap.add_argument("--sample_nmax", type=int, default=40)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    t0 = time.time()
    report = []

    print("part 1: cut convention (q,c) vs (u,v), all pi, 4 <= n <= 6")
    for n in range(4, 7):
        check_cut_convention(n)
        print("  n=%d ok" % n)

    print("part 1b: corollaries 3 (discrepancy) and 4 (mixed difference), all pi, 4 <= n <= 6")
    for n in range(4, 7):
        check_corollaries(n)
        print("  n=%d ok" % n)

    print("part 2: exhaustive ladder (all pi with pi(0)=0 = all toric classes)")
    for n in range(4, args.nmax + 1):
        bound_I = ((n - 1) ** 2) // 4
        maxI = -1
        maxRA = -1
        eqI = 0
        eqRA = 0
        maxTot = -1
        cnt = 0
        full = n <= 8
        for tail in itertools.permutations(range(1, n)):
            pi = (0,) + tail
            I, RA, tot = check_perm(pi, full=full)
            cnt += 1
            maxI = max(maxI, I)
            maxRA = max(maxRA, RA)
            maxTot = max(maxTot, tot)
            if I == bound_I:
                eqI += 1
            if RA == n * bound_I:
                eqRA += 1
        assert maxI <= bound_I, ("H13-I violated", n)
        assert maxRA <= n * bound_I, ("RA violated", n)
        refl_tot = n * n * (n - 1) * (2 * n - 1) // 6
        assert maxTot == refl_tot, ("max Ibar is not the reflection value", n, maxTot, refl_tot)
        line_s = ("  n=%2d classes=%6d  maxI=%d (bound %d, attained %d times)  "
                  "maxRA=%d (n*bound %d, attained %d)  max n^2*Ibar=%d = reflection value"
                  % (n, cnt, maxI, bound_I, eqI, maxRA, n * bound_I, eqRA, maxTot))
        print(line_s)
        report.append(line_s.strip())

    print("part 3: samples up to n=%d (lemmas 1-5 and the theorem)" % args.sample_nmax)
    for n in range(args.nmax + 1, args.sample_nmax + 1):
        bound_I = ((n - 1) ** 2) // 4
        worst = -1
        for pi in samples(n, args.seed):
            I, RA, tot = check_perm(pi, full=False)
            assert I <= bound_I, ("H13-I violated on a sample", n, pi, I)
            assert RA <= n * bound_I, ("RA violated on a sample", n, pi, RA)
            worst = max(worst, I)
        if n % 8 == 0 or n == args.sample_nmax:
            print("  n=%2d ok, max I on samples = %d (bound %d)" % (n, worst, bound_I))
    report.append("samples 9..%d: lemmas hold, I <= floor((n-1)^2/4) on all sampled pi"
                  % args.sample_nmax)

    print("PASS  (%.1f s)" % (time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
