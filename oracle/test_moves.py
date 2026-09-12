"""Mandatory tests for oracle/moves.py (PLAN.md §9.1)."""

import itertools
import random

from moves import (apply_move, apply_word, freely_reduce, delta, identity,
                   sigma, rev, MOVES)


def test_xl_example():
    # PROBLEM.md §3 mandatory example
    assert apply_move((2, 0, 1), "X") == (0, 2, 1)
    assert apply_word((2, 0, 1), "XL") == (2, 1, 0)


def test_problem_word_example():
    # PROBLEM.md §2: XRRXRR from id_4 to sigma_4
    assert apply_word(identity(4), "XRRXRR") == (1, 0, 3, 2) == sigma(4)


def test_inverses_and_orders():
    for n in range(2, 8):
        for p in itertools.permutations(range(n)):
            assert apply_word(p, "LR") == p
            assert apply_word(p, "RL") == p
            assert apply_word(p, "XX") == p
            assert apply_word(p, "L" * n) == p
            for m in MOVES:
                q = apply_move(p, m)
                assert sorted(q) == list(range(n))
            if n > 5:
                break  # exhaustive up to n=5, spot-check beyond


def test_freely_reduce():
    rng = random.Random(12345)
    for n in (4, 5, 6):
        for _ in range(200):
            w = "".join(rng.choice(MOVES) for _ in range(rng.randrange(0, 30)))
            r = freely_reduce(w)
            for p in (identity(n), sigma(n), rev(n)):
                assert apply_word(p, w) == apply_word(p, r)
            assert "XX" not in r and "LR" not in r and "RL" not in r


def test_delta():
    assert delta(0, 3, 6) == 3
    assert delta(0, 4, 6) == 2
    assert delta(5, 0, 6) == 1
    assert delta(2, 2, 7) == 0


def test_sigma_shape():
    assert sigma(4) == (1, 0, 3, 2)
    assert sigma(5) == (1, 0, 4, 3, 2)
    assert sigma(8) == (1, 0, 7, 6, 5, 4, 3, 2)
    assert rev(4) == (3, 2, 1, 0)
    # sigma is an involution
    for n in range(4, 10):
        s = sigma(n)
        assert tuple(s[s[i]] for i in range(n)) == identity(n)


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_"):
            fn()
            print(f"PASS {name}")
    print("ALL TESTS PASS")
