#!/usr/bin/env python3
"""check_H13I.py — independent checker for the session-9 results on H13-I.

H13-I: I(pi) = min over all n^2 double cuts (a, b) of inv(a, b) is at most
floor((n-1)^2/4) for every permutation pi of Z_n, n >= 4.

Everything here is recomputed from scratch (plain O(n^2) inversion counting,
no incremental updates), so the checker is independent of
experiments/line_class_scan.c and experiments/line_profile.c.

Parts:
  1  equivalent numeric forms of the bound (all 4 <= n <= 500);
  2  toric invariance of I and the class reduction to pi(0) = 0 (4 <= n <= 7);
  3  the O(1) update identities used by the C scanner (4 <= n <= 6, all cuts);
  4  max I = floor((n-1)^2/4), attained exactly on the reflection class
     (4 <= n <= 8, exhaustive);
  5  arc criterion (S) is sufficient, and it holds for every extremal pi
     (4 <= n <= 7, exhaustive);
  6  refutation of the families F(m, t) ("line position m carries value t")
     (4 <= n <= 6, exhaustive);
  7  refutation of line-averaging over the cut torus (n = 5, pi = (0,2,4,1,3));
  8  refutation of the Diaconis-Graham route min (D - T) (n = 8, pi(i) = 3i);
  9  refutation of four-direction local minimality (n = 11, explicit cut);
 10  tie to the certified distance tables: d(pi) <= I(pi) + floor(n^2/4)
     (C35) for 4 <= n <= 8 from data/tables/dist_n{n}.bin;
 11  averaging over b at a fixed a is useless on every reflection at every a
     (4 <= n <= 12).

Usage: python3 checks/check_H13I.py [--nmax-full 8]
Version check_H13I-1.0.  Core: no moves needed (pure combinatorics of lines);
part 10 uses the tables certified in session 1 (C2v).
"""
import argparse
import itertools
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# ---------------------------------------------------------------- primitives

def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def line(pi, n, a, b):
    """Line of the double cut (a, b): position x_i = (i-a) % n, value (pi(i)-b) % n."""
    w = [0] * n
    for i in range(n):
        w[(i - a) % n] = (pi[i] - b) % n
    return w


def inv_cut(pi, n, a, b):
    return inversions(line(pi, n, a, b))


def I_of(pi, n):
    return min(inv_cut(pi, n, a, b) for a in range(n) for b in range(n))


def bound(n):
    return (n - 1) ** 2 // 4


def reps(n):
    """Toric class representatives: all pi with pi(0) = 0."""
    for rest in itertools.permutations(range(1, n)):
        yield (0,) + rest


def fail(part, msg):
    print(f"FAIL [{part}]: {msg}")
    sys.exit(1)


# ------------------------------------------------------------------- part 1

def part1():
    for n in range(4, 501):
        c2 = n * (n - 1) // 2
        if c2 - bound(n) != n * n // 4:
            fail(1, f"n={n}: C(n,2) - floor((n-1)^2/4) != floor(n^2/4)")
        if c2 - 2 * bound(n) != n // 2:
            fail(1, f"n={n}: C(n,2) - 2*floor((n-1)^2/4) != floor(n/2)")
    print("part 1 OK: inv <= floor((n-1)^2/4)  <=>  noninv >= floor(n^2/4)"
          "  <=>  Kendall K >= floor(n/2), for 4 <= n <= 500")


# ------------------------------------------------------------------- part 2

def part2(nmax=7):
    for n in range(4, nmax + 1):
        # I is constant on the toric class, and every class meets {pi(0) = 0}
        seen = {}
        for pi in itertools.permutations(range(n)):
            seen[pi] = I_of(pi, n)
        for pi in itertools.permutations(range(n)):
            for al in range(n):
                for be in range(n):
                    pi2 = tuple((pi[(i - al) % n] + be) % n for i in range(n))
                    if seen[pi2] != seen[pi]:
                        fail(2, f"n={n}: I not toric-invariant at {pi}, ({al},{be})")
        gmax = max(seen.values())
        rmax = max(seen[pi] for pi in reps(n))
        if gmax != rmax:
            fail(2, f"n={n}: max over classes {rmax} != max over all pi {gmax}")
    print(f"part 2 OK: I toric-invariant and max over pi(0)=0 = max over all pi, "
          f"4 <= n <= {nmax}")


# ------------------------------------------------------------------- part 3

def part3(nmax=6):
    for n in range(4, nmax + 1):
        for pi in itertools.permutations(range(n)):
            pinv = [0] * n
            for i, v in enumerate(pi):
                pinv[v] = i
            for a in range(n):
                for b in range(n):
                    v = inv_cut(pi, n, a, b)
                    da = inv_cut(pi, n, (a + 1) % n, b) - v
                    db = inv_cut(pi, n, a, (b + 1) % n) - v
                    if da != (n - 1) - 2 * ((pi[a] - b) % n):
                        fail(3, f"n={n} pi={pi} ({a},{b}): a-shift identity")
                    if db != (n - 1) - 2 * ((pinv[b] - a) % n):
                        fail(3, f"n={n} pi={pi} ({a},{b}): b-shift identity")
                    # joint (diagonal) shift
                    dj = inv_cut(pi, n, (a + 1) % n, (b + 1) % n) - v
                    w = line(pi, n, a, b)
                    al, be = w[0], w.index(0)
                    want = 0 if al == 0 else 2 * (n - al - be)
                    if dj != want:
                        fail(3, f"n={n} pi={pi} ({a},{b}): diagonal identity "
                                f"{dj} != {want}")
    print(f"part 3 OK: a-shift, b-shift and joint-shift identities, "
          f"4 <= n <= {nmax}, all pi, all cuts")


# ------------------------------------------------------------------- part 4

def part4(nmax=8):
    for n in range(4, nmax + 1):
        best, args = -1, []
        for pi in reps(n):
            v = I_of(pi, n)
            if v > best:
                best, args = v, [pi]
            elif v == best:
                args.append(pi)
        if best != bound(n):
            fail(4, f"n={n}: max I = {best} != {bound(n)}")
        refl = tuple((-i) % n for i in range(n))
        if args != [refl]:
            fail(4, f"n={n}: extremal reps {args} != [{refl}]")
    print(f"part 4 OK: max I = floor((n-1)^2/4), attained exactly on the "
          f"reflection class, 4 <= n <= {nmax}")


# ------------------------------------------------------------------- part 5

def has_arc(pi, n, k):
    """Does pi map some position arc of length k onto some value arc of length k?"""
    for a in range(n):
        s = set(pi[(a + t) % n] for t in range(k))
        for b in range(n):
            if all(((b + t) % n) in s for t in range(k)):
                return True
    return False


def part5(nmax=7):
    for n in range(4, nmax + 1):
        k = n // 2
        for pi in reps(n):
            v = I_of(pi, n)
            if has_arc(pi, n, k) and v > bound(n):
                fail(5, f"n={n} pi={pi}: (S) holds but I = {v} > {bound(n)}")
            if v == bound(n) and not has_arc(pi, n, k):
                fail(5, f"n={n} pi={pi}: extremal but (S) fails")
    print(f"part 5 OK: arc criterion (S) sufficient and satisfied by every "
          f"extremal pi, 4 <= n <= {nmax}")


# ------------------------------------------------------------------- part 6

def part6(nmax=6):
    for n in range(4, nmax + 1):
        worst = {}
        for m in range(n):
            for t in range(n):
                worst[(m, t)] = -1
        for pi in reps(n):
            for (m, t) in worst:
                v = min(inv_cut(pi, n, a, (pi[(a + m) % n] - t) % n)
                        for a in range(n))
                if v > worst[(m, t)]:
                    worst[(m, t)] = v
        good = [f for f in worst if worst[f] <= bound(n)]
        if good:
            fail(6, f"n={n}: families {good} would suffice after all")
    print(f"part 6 OK: no family F(m, t) of n cuts suffices, 4 <= n <= {nmax}")


# ------------------------------------------------------------------- part 7

def min_line_sum(pi, n):
    """Smallest sum of inv over a line of the cut torus, over all n+1 directions:
    (1, g) for g = 0 .. n-1 and (0, 1)."""
    best = min(sum(inv_cut(pi, n, a, (a * g + t) % n) for a in range(n))
               for g in range(n) for t in range(n))
    best = min(best, min(sum(inv_cut(pi, n, t, b) for b in range(n))
                         for t in range(n)))
    return best


def part7():
    n, pi = 5, (0, 2, 1, 4, 3)
    best = min_line_sum(pi, n)
    if best <= n * bound(n):
        fail(7, f"line averaging gives {best}/{n} <= {bound(n)} after all")
    # exhaustive at n = 5: every pi is covered and this one is the worst
    worst = max(min_line_sum(p, n) for p in reps(n))
    if worst != best:
        fail(7, f"n=5: worst line-average {worst} != {best}")
    print(f"part 7 OK: averaging over any single line of the cut torus fails "
          f"already at n=5, worst pi={pi}: best line average {best}/{n} > "
          f"{bound(n)} (exhaustive over all classes at n=5)")


# ------------------------------------------------------------------- part 8

def cayley(w):
    n = len(w)
    seen, c = [False] * n, 0
    for i in range(n):
        if not seen[i]:
            c += 1
            j = i
            while not seen[j]:
                seen[j] = True
                j = w[j]
    return n - c


def part8():
    n, pi = 8, tuple((3 * i) % 8 for i in range(8))
    best = None
    for a in range(n):
        for b in range(n):
            w = line(pi, n, a, b)
            v = sum(abs(j - w[j]) for j in range(n)) - cayley(w)
            if best is None or v < best:
                best = v
    if best <= bound(n):
        fail(8, f"D - T route gives {best} <= {bound(n)} after all")
    # and the route is sound where it applies: inv <= D - T at every cut
    for a in range(n):
        for b in range(n):
            w = line(pi, n, a, b)
            if inversions(w) > sum(abs(j - w[j]) for j in range(n)) - cayley(w):
                fail(8, "Diaconis-Graham inequality violated (implementation bug)")
    print(f"part 8 OK: min (D - T) = {best} > {bound(n)} at n=8, pi(i)=3i — "
          f"the Diaconis-Graham route cannot prove H13-I")


# ------------------------------------------------------------------- part 9

def part9():
    n, pi = 11, (0, 2, 10, 6, 5, 1, 9, 8, 7, 4, 3)
    a, b = 4, 1
    v = inv_cut(pi, n, a, b)
    if v != 26 or bound(n) != 25:
        fail(9, f"n=11 counterexample: inv={v}, bound={bound(n)}")
    for d in ((1, 0), (0, 1), (1, 1), (1, -1)):
        for k in range(1, n):
            if inv_cut(pi, n, (a + k * d[0]) % n, (b + k * d[1]) % n) < v:
                fail(9, f"cut ({a},{b}) is not locally minimal along {d}")
    if I_of(pi, n) != 18:
        fail(9, "I(pi) != 18 for the n=11 counterexample")
    print(f"part 9 OK: n=11, pi={pi}, cut ({a},{b}): inv = 26 > 25 and the cut "
          f"is locally minimal along (1,0), (0,1), (1,1), (1,-1) — four-direction "
          f"local minimality does not imply H13-I (I(pi) = 18)")


# ------------------------------------------------------------------ part 10

def rank_perm(p):
    n = len(p)
    fact = [1] * (n + 1)
    for i in range(1, n + 1):
        fact[i] = fact[i - 1] * i
    r = 0
    for i in range(n):
        smaller = sum(1 for j in range(i + 1, n) if p[j] < p[i])
        r += smaller * fact[n - 1 - i]
    return r


def part10(nmax=8):
    done = []
    for n in range(4, nmax + 1):
        path = os.path.join(ROOT, "data", "tables", f"dist_n{n}.bin")
        if not os.path.exists(path):
            continue
        with open(path, "rb") as f:
            dist = f.read()
        for pi in itertools.permutations(range(n)):
            d = dist[rank_perm(list(pi))]
            v = I_of(pi, n)
            if d > v + n * n // 4:
                fail(10, f"n={n} pi={pi}: d={d} > I+floor(n^2/4)={v + n*n//4}")
            if v > bound(n):
                fail(10, f"n={n} pi={pi}: I={v} > {bound(n)}")
        done.append(n)
    print(f"part 10 OK: d(pi) <= I(pi) + floor(n^2/4) and I(pi) <= "
          f"floor((n-1)^2/4) against the certified tables, n in {done}")


# ------------------------------------------------------------------ part 11

def part11(nmax=12):
    """Averaging over b at a fixed position cut a is useless on every reflection
    at every a: Q(a) = sum over pairs in the a-order of (pi(q)-pi(p)) mod n is
    constant and the resulting bound noninv >= C(n,2) - Q(a)/n equals
    (n^2-1)/6 < floor(n^2/4)."""
    for n in range(4, nmax + 1):
        pred = n * (n * (n - 1) // 2) - n * (n * n - 1) // 6
        for h in range(n):
            pi = [(h - i) % n for i in range(n)]
            for a in range(n):
                order = sorted(range(n), key=lambda i: (i - a) % n)
                q = sum((pi[y] - pi[x]) % n
                        for k, x in enumerate(order) for y in order[k + 1:])
                if q != pred:
                    fail(11, f"n={n} h={h} a={a}: Q={q} != {pred}")
        # the bound it yields, times 6, versus floor(n^2/4) times 6
        if (n * n - 1) >= 6 * (n * n // 4):
            fail(11, f"n={n}: (n^2-1)/6 >= floor(n^2/4), claim as stated is wrong")
    print(f"part 11 OK: on every reflection Q(a) is constant and the "
          f"average-over-b bound gives only noninv >= (n^2-1)/6 < floor(n^2/4), "
          f"4 <= n <= {nmax}")



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax-full", type=int, default=8,
                    help="largest n for the exhaustive part 4")
    ap.add_argument("--skip", default="", help="comma separated part numbers to skip")
    args = ap.parse_args()
    skip = {int(x) for x in args.skip.split(",") if x.strip()}

    parts = [
        (1, part1, ()),
        (2, part2, (7,)),
        (3, part3, (6,)),
        (4, part4, (args.nmax_full,)),
        (5, part5, (7,)),
        (6, part6, (6,)),
        (7, part7, ()),
        (8, part8, ()),
        (9, part9, ()),
        (10, part10, (8,)),
        (11, part11, (12,)),
    ]
    for num, fn, a in parts:
        if num in skip:
            print(f"part {num} skipped")
            continue
        fn(*a)
    print("PASS")


if __name__ == "__main__":
    main()
