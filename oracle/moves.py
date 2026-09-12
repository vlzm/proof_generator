"""Reference implementation of the LRX moves.

Conventions (PROBLEM.md §3):
- states are tuples of ints 0..n-1, 0-indexed;
- words are executed left to right: "XL" means first X, then L;
- L: (x0, x1, ..., x_{n-1}) -> (x1, ..., x_{n-1}, x0);
- R: inverse of L;
- X: swap the first two entries.

Mandatory check: apply_word((2,0,1), "XL") == (2,1,0).
"""

CORE_VERSION = "oracle-1.0"

MOVES = ("L", "R", "X")

INVERSE = {"L": "R", "R": "L", "X": "X"}


def apply_move(state, move):
    """Apply a single move to a state (tuple), returning a new tuple."""
    if move == "L":
        return state[1:] + state[:1]
    if move == "R":
        return state[-1:] + state[:-1]
    if move == "X":
        return (state[1], state[0]) + state[2:]
    raise ValueError(f"unknown move: {move!r}")


def apply_word(state, word):
    """Apply a word (string over L, R, X) left to right."""
    for m in word:
        state = apply_move(state, m)
    return state


def freely_reduce(word):
    """Cancel adjacent XX, LR, RL pairs until none remain."""
    out = []
    for m in word:
        if m not in MOVES:
            raise ValueError(f"unknown move: {m!r}")
        if out and INVERSE[out[-1]] == m:
            out.pop()
        else:
            out.append(m)
    return "".join(out)


def delta(a, b, n):
    """Circular distance between positions a and b on Z_n."""
    d = (b - a) % n
    return min(d, n - d)


def identity(n):
    return tuple(range(n))


def sigma(n):
    """The conjectured unique farthest vertex sigma_n[i] = (1-i) mod n."""
    return tuple((1 - i) % n for i in range(n))


def rev(n):
    """Plain reversal (n-1, n-2, ..., 0)."""
    return tuple(range(n - 1, -1, -1))


def word_stats(word):
    """Return (length, N_X, N_rot) for a word."""
    n_x = word.count("X")
    return len(word), n_x, len(word) - n_x
