"""h13_pair_corner-1.0

Session 9: exploration around H13-I (docs/notes/h13_line_model.md §6).

Part 1 (--verify-formula): cross-checks the pair-corner formula
f(d,e) = n(d+e) - 2*d*e  (proof: docs/proofs/C37_pair_corner_formula.md)
against a direct brute-force count over all n^2 corners, for random
permutations.

Part 2 (--scan): extends the H13-I bound check (I(pi) <= floor((n-1)^2/4))
beyond the exhaustive range (4 <= n <= 10, see C33) to structured families
and random samples at larger n. This is error-hunting (AGENTS.md rule 9),
not a proof of the general case.

Part 3 (--hillclimb): local search (random-restart hill climbing on
transpositions) trying to find a permutation violating the H13-I bound,
at n where exhaustive search is not affordable (n >= 11).

Usage:
    python3 experiments/h13_pair_corner.py --verify-formula
    python3 experiments/h13_pair_corner.py --scan --sizes 15,20,30,50
    python3 experiments/h13_pair_corner.py --hillclimb --sizes 11,12,13,14 \
        --iters 1500 --restarts 5
"""
import argparse
import random
import time


def inv_count(w):
    n = len(w)
    if n <= 1:
        return 0

    def sort_count(a):
        if len(a) <= 1:
            return a, 0
        mid = len(a) // 2
        left, cl = sort_count(a[:mid])
        right, cr = sort_count(a[mid:])
        merged = []
        i = j = 0
        inv = cl + cr
        while i < len(left) and j < len(right):
            if left[i] <= right[j]:
                merged.append(left[i])
                i += 1
            else:
                merged.append(right[j])
                j += 1
                inv += len(left) - i
        merged.extend(left[i:])
        merged.extend(right[j:])
        return merged, inv

    _, inv = sort_count(w)
    return inv


def I_of_pi(pi):
    """min over all n^2 corners (a,b) of inv(w); see PROBLEM/h13_line_model.md."""
    n = len(pi)
    best = None
    for a in range(n):
        rotated = [pi[(a + j) % n] for j in range(n)]
        for b in range(n):
            w = [(x - b) % n for x in rotated]
            v = inv_count(w)
            if best is None or v < best:
                best = v
            if best == 0:
                return 0
    return best


def f_formula(d, e, n):
    return n * (d + e) - 2 * d * e


def discordant_count_direct(pi, i1, i2, n):
    y1, y2 = pi[i1], pi[i2]
    cnt = 0
    for a in range(n):
        for b in range(n):
            j1 = (i1 - a) % n
            j2 = (i2 - a) % n
            w1 = (y1 - b) % n
            w2 = (y2 - b) % n
            if (j1 < j2) != (w1 < w2):
                cnt += 1
    return cnt


def verify_formula(sizes, seed=1):
    rng = random.Random(seed)
    total_pairs = 0
    mismatches = 0
    for n in sizes:
        pi = list(range(n))
        rng.shuffle(pi)
        for i1 in range(n):
            for i2 in range(i1 + 1, n):
                d = i2 - i1
                e = (pi[i2] - pi[i1]) % n
                fd = f_formula(d, e, n)
                direct = discordant_count_direct(pi, i1, i2, n)
                total_pairs += 1
                if fd != direct:
                    mismatches += 1
                    print(f"MISMATCH n={n} i1={i1} i2={i2} d={d} e={e} "
                          f"formula={fd} direct={direct}")
        print(f"n={n}: pi={pi} — {n*(n-1)//2} pairs checked")
    print(f"TOTAL pairs checked: {total_pairs}, mismatches: {mismatches}")
    return mismatches == 0


def scan(sizes, seed=42):
    rng = random.Random(seed)
    all_ok = True
    for n in sizes:
        bound = ((n - 1) ** 2) // 4
        families = []
        families.append((f"reflection sigma_n", [(-i) % n for i in range(n)]))
        families.append((f"reverse", list(range(n))[::-1]))
        for a, b in [(2, 1), (3, 5), (n // 2 + 1, 3)]:
            pi = [(a * i + b) % n for i in range(n)]
            if len(set(pi)) == n:
                families.append((f"affine a={a} b={b}", pi))
        for t in range(5):
            pi = list(range(n))
            rng.shuffle(pi)
            families.append((f"random#{t}", pi))
        for label, pi in families:
            I = I_of_pi(pi)
            ok = I <= bound
            all_ok &= ok
            status = "OK" if ok else "VIOLATION"
            print(f"n={n:4d} I={I:6d} bound={bound:6d} {status}  {label}")
    return all_ok


def hillclimb(sizes, iters, restarts, seed=0):
    all_ok = True
    for n in sizes:
        rng = random.Random(seed + n)
        bound = ((n - 1) ** 2) // 4
        global_best = -1
        global_best_pi = None
        for r in range(restarts):
            pi = list(range(n))
            rng.shuffle(pi)
            cur = I_of_pi(pi)
            for _ in range(iters):
                i, j = rng.sample(range(n), 2)
                pi[i], pi[j] = pi[j], pi[i]
                val = I_of_pi(pi)
                if val >= cur:
                    cur = val
                else:
                    pi[i], pi[j] = pi[j], pi[i]
            if cur > global_best:
                global_best = cur
                global_best_pi = pi[:]
        ok = global_best <= bound
        all_ok &= ok
        status = "VIOLATION!" if not ok else ("tight" if global_best == bound else "ok")
        print(f"n={n} bound={bound} hillclimb_max_I={global_best} {status} "
              f"pi={global_best_pi}")
    return all_ok


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--verify-formula", action="store_true")
    p.add_argument("--scan", action="store_true")
    p.add_argument("--hillclimb", action="store_true")
    p.add_argument("--sizes", default="15,20,30,50")
    p.add_argument("--iters", type=int, default=1500)
    p.add_argument("--restarts", type=int, default=5)
    p.add_argument("--seed", type=int, default=1)
    args = p.parse_args()

    sizes = [int(x) for x in args.sizes.split(",") if x]
    t0 = time.time()
    ok = True
    if args.verify_formula:
        ok &= verify_formula([5, 6, 7], seed=args.seed)
    if args.scan:
        ok &= scan(sizes, seed=args.seed)
    if args.hillclimb:
        ok &= hillclimb(sizes, args.iters, args.restarts, seed=args.seed)
    print(f"elapsed: {time.time()-t0:.1f}s  overall: {'OK' if ok else 'VIOLATION FOUND'}")


if __name__ == "__main__":
    main()
