"""delta_lemma_check-1.0

Independent check of the increment lemma used in the session-9 attempt on
H13-I (see docs/notes/h13_line_model.md, section "Session 9"):

    inv_{c+1}(sigma) - inv_c(sigma) = n - 1 - 2*k_c

where sigma is a permutation of {0,...,n-1} (values), inv_c(sigma) is the
number of inversions of the sequence (sigma_i - c) mod n read left to
right, and k_c is the position (0-indexed) of value c in sigma
(k_c = sigma^{-1}(c)).

This is a pure statement about permutations and inversions, unrelated to
the LRX oracle, so it is checked directly against a brute-force inversion
count rather than against oracle/moves.py.

Ladder (AGENTS.md rule 9): n = 4..9 exhaustive; n = 10, 20, 50, 100, 101
random samples with fixed seed.
"""
import itertools
import random
import sys


def inv_c(sigma, c, n):
    w = [(x - c) % n for x in sigma]
    inv = 0
    for i in range(n):
        wi = w[i]
        for j in range(i + 1, n):
            if wi > w[j]:
                inv += 1
    return inv


def check_one(sigma, n):
    pos_of = [0] * n
    for idx, v in enumerate(sigma):
        pos_of[v] = idx
    base = inv_c(sigma, 0, n)
    cur = base
    for c in range(n):
        nxt = inv_c(sigma, (c + 1) % n, n)
        actual = nxt - cur
        predicted = n - 1 - 2 * pos_of[c]
        if actual != predicted:
            return False, (sigma, c, actual, predicted)
        cur = nxt
    # closed walk: inv_n == inv_0
    if cur != base:
        return False, (sigma, "closed", cur, base)
    return True, None


def exhaustive(n):
    checked = 0
    for perm in itertools.permutations(range(n)):
        ok, info = check_one(list(perm), n)
        checked += 1
        if not ok:
            return checked, info
    return checked, None


def sampled(n, trials, seed):
    rng = random.Random(seed)
    checked = 0
    base = list(range(n))
    for _ in range(trials):
        rng.shuffle(base)
        ok, info = check_one(list(base), n)
        checked += 1
        if not ok:
            return checked, info
    return checked, None


def main():
    print("n, mode, checked, result")
    for n in range(4, 10):
        checked, info = exhaustive(n)
        print(f"{n}, exhaustive, {checked}, {'PASS' if info is None else info}")
        if info is not None:
            sys.exit(1)
    for n in (10, 20, 50, 100, 101):
        trials = 20000 if n <= 20 else 2000
        checked, info = sampled(n, trials, seed=1)
        print(f"{n}, sampled(seed=1), {checked}, {'PASS' if info is None else info}")
        if info is not None:
            sys.exit(1)
    print("ALL PASS")


if __name__ == "__main__":
    main()
