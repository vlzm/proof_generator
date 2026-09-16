"""checks/check_C37_C38.py — independent Python re-implementation checking
claims C37 (Lemma D: exact O(n^2) reduction formula for I(pi) via n position
rotations) and C38 (Lemma E: the averaged bound avg_a h(a) <= floor((n-1)^2/4)),
plus the closed-form identity for T1 = sum_a inv(beta_a) used in
docs/notes/h13_line_model.md §7.

Independent from experiments/h13i_avg_bound.c: written from scratch in Python,
against the *definition* of I(pi) via all n^2 double cuts (matching
experiments/line_model.py's convention), not against the C code.

Usage: python3 checks/check_C37_C38.py [nmax]   (default nmax=8; exhaustive
over all n! permutations for 3 <= n <= nmax; keep nmax <= 9 for a check to
finish in reasonable time in pure Python).

Version check_C37_C38-1.0.
"""
import itertools
import sys


def inv_seq(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def brute_I(pi):
    """I(pi) = min over all n^2 double cuts (q, c) of inv of the relabelled
    line, exactly as in experiments/line_model.py."""
    n = len(pi)
    best = None
    for q in range(n):
        for c in range(n):
            shift = (q + 1 - c) % n
            w = [(pi[(q + 1 + j) % n] - shift) % n for j in range(n)]
            iv = inv_seq(w)
            if best is None or iv < best:
                best = iv
    return best


def h_of_rotation(pi, a):
    """h(a) = min_b inv(shift_b(beta_a)) via Lemma D's closed form."""
    n = len(pi)
    beta = [pi[(k + a + 1) % n] for k in range(n)]
    pos = [0] * n
    for k, v in enumerate(beta):
        pos[v] = k
    inv_beta = inv_seq(beta)
    P = 0
    best_D = 0  # k = 0
    for k in range(1, n + 1):
        P += pos[k - 1]
        D = P - k * (n - 1) / 2
        if D > best_D:
            best_D = D
    return inv_beta - 2 * best_D


def lemma_d_I(pi):
    n = len(pi)
    return min(h_of_rotation(pi, a) for a in range(n))


def T1_direct(pi):
    n = len(pi)
    return sum(inv_seq([pi[(k + a + 1) % n] for k in range(n)]) for a in range(n))


def T1_formula(pi):
    n = len(pi)
    inv0 = inv_seq(pi)
    W = sum(j - i for i in range(n) for j in range(i + 1, n) if pi[i] > pi[j])
    return n * inv0 - 2 * W + n * (n * n - 1) // 6


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    lemma_d_fail = 0
    lemma_e_fail = 0
    t1_fail = 0
    for n in range(3, nmax + 1):
        bound = (n - 1) ** 2 // 4
        max_avg = -1.0
        argmax = None
        cnt = 0
        for pi in itertools.permutations(range(n)):
            pi = list(pi)
            bi = brute_I(pi)
            di = lemma_d_I(pi)
            if abs(bi - di) > 1e-9:
                print(f"C37 (Lemma D) MISMATCH n={n} pi={pi} brute={bi} lemma_d={di}")
                lemma_d_fail += 1
            if T1_direct(pi) != T1_formula(pi):
                print(f"T1 formula MISMATCH n={n} pi={pi}")
                t1_fail += 1
            avg_h = sum(h_of_rotation(pi, a) for a in range(n)) / n
            if avg_h > max_avg + 1e-9:
                max_avg = avg_h
                argmax = pi
                cnt = 1
            elif abs(avg_h - max_avg) <= 1e-9:
                cnt += 1
        exceeds = max_avg > bound + 1e-9
        if exceeds:
            lemma_e_fail += 1
        print(f"n={n}: C37 mismatches={lemma_d_fail}, T1 mismatches={t1_fail}, "
              f"C38 max avg_a h(a)={max_avg:.4f} vs bound {bound} "
              f"(exceeds={exceeds}, #argmax={cnt}, one argmax={argmax})")
    status = "PASS" if (lemma_d_fail == 0 and lemma_e_fail == 0 and t1_fail == 0) else "FAIL"
    print(status)
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
