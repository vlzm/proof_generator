"""H13-I exploration: existence of a double cut (q, c) with few inversions.

Definitions match docs/notes/h13_line_model.md §0. For pi in S_n and a cut
(q, c), the line w_j = (pi(q+1+j) - (q+1-c)) mod n, j = 0..n-1, has inv(w)
inversions. I(pi) = min over the n^2 cuts of inv(w). H13-I conjectures
I(pi) <= floor((n-1)^2/4) for all pi, n >= 4.

This script:
1. Checks the closed-form lemma f(d, e) = n(d+e) - 2*d*e: for a fixed
   unordered pair of positions with position-gap d and value-gap e
   (both in {1,...,n-1}, taken from either ordering), the pair is inverted
   at exactly f(d, e) of the n^2 cuts. (PROVED; see report and CLAIMS C37.)
2. Exhaustively verifies max_pi I(pi) = floor((n-1)^2/4), attained on
   exactly n permutations (the reflections pi_h(i) = (h - i) mod n).
3. Checks three restricted (single-degree-of-freedom) families that were
   candidates for an easy existence proof, all REFUTED as insufficient:
   - rotation-only (c fixed at a canonical value, only q varies),
   - "front pivot" diagonal (for each q, c chosen so w_0 = 0),
   - best single row (choose q minimizing the row total S(q), then use
     the row average bound) -- fails already on the identity permutation.
"""
import itertools


def inv_count(seq):
    n = len(seq)
    c = 0
    for i in range(n):
        for j in range(i + 1, n):
            if seq[i] > seq[j]:
                c += 1
    return c


def f_formula(n, d, e):
    return n * (d + e) - 2 * d * e


def all_cuts_inv(pi):
    n = len(pi)
    table = {}
    for q in range(n):
        pos_order = [(q + 1 + j) % n for j in range(n)]
        v = [pi[p] for p in pos_order]
        for c in range(n):
            t = (q + 1 - c) % n
            w = [(x - t) % n for x in v]
            table[(q, c)] = inv_count(w)
    return table


def I_of_pi(pi):
    return min(all_cuts_inv(pi).values())


def rotation_only_min(pi):
    """c fixed so that t = 0 (no value shift); only q varies (n candidates)."""
    n = len(pi)
    best = None
    for q in range(n):
        seq = [pi[(q + 1 + j) % n] for j in range(n)]
        iv = inv_count(seq)
        if best is None or iv < best:
            best = iv
    return best


def front_pivot_min(pi):
    """For each q, choose c so w_0 = 0 (n candidates, the 'diagonal' family)."""
    n = len(pi)
    best = None
    for q in range(n):
        pos_order = [(q + 1 + j) % n for j in range(n)]
        v = [pi[p] for p in pos_order]
        t = v[0]
        w = [(x - t) % n for x in v]
        iv = inv_count(w)
        if best is None or iv < best:
            best = iv
    return best


def row_sum(pi, q):
    n = len(pi)
    pos_order = [(q + 1 + j) % n for j in range(n)]
    v = [pi[p] for p in pos_order]
    total = 0
    for i in range(n):
        for j in range(i + 1, n):
            total += (v[j] - v[i]) % n
    return total


def best_row_average_bound(pi):
    """min_q S(q) / n, floored -- the pigeonhole bound from the best single row."""
    n = len(pi)
    m = min(row_sum(pi, q) for q in range(n))
    return m // n


def verify_f_formula(n, trials=200, seed=0):
    import random
    rng = random.Random(seed)
    for _ in range(trials):
        pi = list(range(n))
        rng.shuffle(pi)
        i, j = rng.sample(range(n), 2)
        d = (j - i) % n
        e = (pi[j] - pi[i]) % n
        cnt = 0
        for q in range(n):
            pos_order = [(q + 1 + k) % n for k in range(n)]
            rank = {p: idx for idx, p in enumerate(pos_order)}
            for c in range(n):
                t = (q + 1 - c) % n
                wi = (pi[i] - t) % n
                wj = (pi[j] - t) % n
                inverted = (wi > wj) if rank[i] < rank[j] else (wj > wi)
                cnt += inverted
        assert cnt == f_formula(n, d, e), (n, i, j, d, e, cnt, f_formula(n, d, e))
    return True


def main():
    print("# f(d,e) formula spot-check")
    for n in range(4, 10):
        verify_f_formula(n, 100, seed=n)
        print(f"n={n}: f(d,e) formula OK (100 random pairs)")

    print()
    print("# exhaustive checks, 4 <= n <= 8")
    header = f"{'n':>2} {'T=floor((n-1)^2/4)':>20} {'max I(pi)':>10} {'#argmax':>8} " \
              f"{'max rot-only':>13} {'max front-pivot':>16} {'max row-avg-bound fails':>24}"
    print(header)
    for n in range(4, 9):
        T = (n - 1) ** 2 // 4
        max_I = -1
        argmax_count = 0
        max_rot = -1
        max_fp = -1
        row_avg_fail = 0
        for pi in itertools.permutations(range(n)):
            table = all_cuts_inv(pi)
            I = min(table.values())
            if I > max_I:
                max_I = I
                argmax_count = 1
            elif I == max_I:
                argmax_count += 1
            r = rotation_only_min(pi)
            if r > max_rot:
                max_rot = r
            fp = front_pivot_min(pi)
            if fp > max_fp:
                max_fp = fp
            if best_row_average_bound(pi) > T:
                row_avg_fail += 1
        print(f"{n:>2} {T:>20} {max_I:>10} {argmax_count:>8} "
              f"{max_rot:>13} {max_fp:>16} {row_avg_fail:>24}")


if __name__ == "__main__":
    main()
