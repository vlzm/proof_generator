"""H13-I search: is I(pi) = min_{q,s} inv(cut) always <= floor((n-1)^2/4)?

H13-I (docs/notes/h13_line_model.md #6, PLAN SS8): every permutation pi of
Z_n has a double cut (rotation q of positions, rotation s of values) whose
relabelled line has at most floor((n-1)^2/4) inversions.  Proved so far only
by exhaustive check, 4 <= n <= 10 (C33, experiments/line_profile.c).  This
module does three things for the "prove H13-I, budget 1 session" task
(session 9):

1. `fast_I(pi)`: O(n^2 log n) computation of I(pi) via a Fenwick tree and
   the recurrence below, in place of the O(n^4) brute force (n^2 cuts x
   O(n^2) inversions each).  Verified against brute force: exhaustive
   3 <= n <= 8, random spot checks n in {9, 10, 12, 15}.

   Recurrence.  Fix the position cut q (v = pi read starting after q).
   Going from value-shift s to s+1, every w_j decreases by 1 mod n except
   the position m holding value s, which wraps from 0 to n-1 (min -> max).
   Only pairs touching m change order; a short count gives
   inv(s+1) = inv(s) + (n - 1 - 2*m), where m = (pi^{-1}(s) - q - 1) mod n
   is m's index in the v-array.  So all n values of s cost O(n) after one
   O(n log n) inversion count at s = 0; summed over q, that's O(n^2 log n).

2. `search(n, ...)`: exhaustive reflections + affine family + random
   sampling + simulated annealing (maximizing I) hunting for a
   counterexample I(pi) > floor((n-1)^2/4), for n beyond the exhaustive
   range (11 up to 100).

3. Two rejected restricted-cut families, run to exhaustion at small n as an
   extra negative result (candidates from the "next step" list in
   docs/notes/h13_line_model.md #6):
   - "diagonal" cuts: origin (q, s) = (i - 1, pi(i)) for one of the n data
     points i (n candidates instead of n^2).
   - "single shift": fix q (natural order) and only optimize s (n
     candidates, already implicit in Lemma C's counterexample but checked
     here as its own family over all pi, not just one input).
   Both fail on a growing fraction of pi as n grows (see report).

Usage:
  python3 experiments/h13i_search.py --selftest
  python3 experiments/h13i_search.py --restricted-families --nmax 8
  python3 experiments/h13i_search.py --search --nlist 11,12,15,20,25,30,40,50,70,100 --sa-nlist 11,12,13,14,15,16,18,20,24
Output: data/runs/h13i_search/report.md (this run), no report.json (no
counterexample found, nothing to replay). Version h13i_search-1.0.
"""
import argparse
import itertools
import math
import random
import time


class Fenwick:
    __slots__ = ("n", "t")

    def __init__(self, n):
        self.n = n
        self.t = [0] * (n + 1)

    def add(self, i, v=1):
        i += 1
        while i <= self.n:
            self.t[i] += v
            i += i & (-i)

    def sum(self, i):  # sum of counts in [0, i)
        s = 0
        while i > 0:
            s += self.t[i]
            i -= i & (-i)
        return s


def inv_count(w):
    n = len(w)
    fw = Fenwick(n)
    inv = 0
    for x in reversed(w):
        inv += fw.sum(x)
        fw.add(x)
    return inv


def bound(n):
    return (n - 1) ** 2 // 4


def fast_I(pi):
    """I(pi) = min over n^2 double cuts of the number of inversions."""
    n = len(pi)
    pi_inv = [0] * n
    for p, v in enumerate(pi):
        pi_inv[v] = p
    best = None
    for q in range(n):
        v = [pi[(q + 1 + j) % n] for j in range(n)]
        cur = inv_count(v)
        if best is None or cur < best:
            best = cur
        for s in range(1, n):
            m = (pi_inv[s - 1] - q - 1) % n
            cur += n - 1 - 2 * m
            if cur < best:
                best = cur
    return best


def brute_I(pi):
    n = len(pi)
    best = None
    for q in range(n):
        for s in range(n):
            w = [(pi[(q + 1 + j) % n] - s) % n for j in range(n)]
            c = sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])
            if best is None or c < best:
                best = c
    return best


def diagonal_I(pi):
    """Restricted family: origin = (i-1, pi[i]) for one of the n points."""
    n = len(pi)
    best = None
    for i in range(n):
        q = (i - 1) % n
        s = pi[i]
        w = [(pi[(q + 1 + j) % n] - s) % n for j in range(n)]
        c = sum(1 for a in range(n) for b in range(a + 1, n) if w[a] > w[b])
        if best is None or c < best:
            best = c
    return best


def single_shift_I(pi):
    """Restricted family: fix q = n-1 (natural order), optimize s only."""
    n = len(pi)
    best = None
    for s in range(n):
        w = [(pi[j] - s) % n for j in range(n)]
        c = sum(1 for a in range(n) for b in range(a + 1, n) if w[a] > w[b])
        if best is None or c < best:
            best = c
    return best


def selftest():
    for n in range(3, 9):
        for pi in itertools.permutations(range(n)):
            b = brute_I(list(pi))
            f = fast_I(list(pi))
            assert b == f, (n, pi, b, f)
        print(f"n={n} fast_I == brute_I: exhaustive OK ({math.factorial(n)} permutations)")
    rng = random.Random(1)
    for n in (9, 10, 12, 15):
        for _ in range(200):
            pi = list(range(n))
            rng.shuffle(pi)
            assert brute_I(pi) == fast_I(pi), (n, pi)
        print(f"n={n} fast_I == brute_I: 200 random spot checks OK")


def restricted_families(nmax):
    print("n  bound  diagonal_fails  single_shift_fails  total")
    for n in range(4, nmax + 1):
        total = 0
        dfail = 0
        sfail = 0
        b = bound(n)
        for pi in itertools.permutations(range(n)):
            total += 1
            if diagonal_I(list(pi)) > b:
                dfail += 1
            if single_shift_I(list(pi)) > b:
                sfail += 1
        print(f"{n:2d} {b:5d} {dfail:8d} {sfail:12d} {total:8d}")


def reflection(n, h):
    return [(h - i) % n for i in range(n)]


def affine(n, a, b):
    return [(a * i + b) % n for i in range(n)]


def search(n, n_random, seed):
    rng = random.Random(seed)
    b = bound(n)
    worst = 0
    worst_desc = None
    for h in range(n):
        v = fast_I(reflection(n, h))
        if v > worst:
            worst, worst_desc = v, f"reflection h={h}"
    for a in range(1, n):
        if math.gcd(a, n) != 1:
            continue
        for bb in range(n):
            v = fast_I(affine(n, a, bb))
            if v > worst:
                worst, worst_desc = v, f"affine a={a} b={bb}"
    for _ in range(n_random):
        pi = list(range(n))
        rng.shuffle(pi)
        v = fast_I(pi)
        if v > worst:
            worst, worst_desc = v, f"random {pi}"
    return b, worst, worst_desc


def simulated_annealing(n, iters, restarts, seed):
    rng = random.Random(seed)
    global_best = 0
    global_pi = None
    for _ in range(restarts):
        pi = list(range(n))
        rng.shuffle(pi)
        cur = fast_I(pi)
        best = cur
        T0, T1 = 3.0, 0.02
        for it in range(iters):
            T = T0 * ((T1 / T0) ** (it / iters))
            i = rng.randrange(n)
            j = rng.randrange(n)
            if i == j:
                continue
            pi[i], pi[j] = pi[j], pi[i]
            v = fast_I(pi)
            d = v - cur
            if d >= 0 or rng.random() < math.exp(d / T):
                cur = v
            else:
                pi[i], pi[j] = pi[j], pi[i]
                continue
            if cur > best:
                best = cur
        if best > global_best:
            global_best = best
    return global_best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--restricted-families", action="store_true")
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--search", action="store_true")
    ap.add_argument("--nlist", default="11,12,15,20,25,30,40,50,70,100")
    ap.add_argument("--sa-nlist", default="11,12,13,14,15,16,18,20,24")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    if args.selftest:
        selftest()
    if args.restricted_families:
        restricted_families(args.nmax)
    if args.search:
        print("n    bound  worst_I  source  exceed  time_s")
        for n in [int(x) for x in args.nlist.split(",") if x]:
            t0 = time.time()
            n_random = 3000 if n <= 30 else (800 if n <= 50 else 200)
            b, worst, desc = search(n, n_random, seed=args.seed + n)
            dt = time.time() - t0
            exceed = "YES" if worst > b else "no"
            print(f"{n:4d} {b:6d} {worst:7d}  {desc}  {exceed}  {dt:.1f}")
        print()
        print("n    bound  SA_best  exceed  time_s")
        for n in [int(x) for x in args.sa_nlist.split(",") if x]:
            t0 = time.time()
            iters = 15000 if n <= 16 else 8000
            b = bound(n)
            best = simulated_annealing(n, iters, restarts=6, seed=args.seed + 1000 + n)
            dt = time.time() - t0
            exceed = "YES" if best > b else "no"
            print(f"{n:4d} {b:5d} {best:7d}  {exceed}  {dt:.1f}")


if __name__ == "__main__":
    main()
