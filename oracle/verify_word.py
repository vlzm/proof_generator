"""Word checker (PLAN.md §9.1): validates alphabet and the result of a word.

verify_word(start, word, expected) -> (ok, info)
  - checks word alphabet is a subset of {L, R, X};
  - applies the word left to right with reference moves;
  - reports length, N_X, N_rot, and the freely reduced length.
"""

from moves import apply_word, freely_reduce, word_stats, MOVES


def verify_word(start, word, expected=None):
    for ch in word:
        if ch not in MOVES:
            return False, {"error": f"illegal letter {ch!r}"}
    result = apply_word(start, word)
    length, n_x, n_rot = word_stats(word)
    info = {
        "result": result,
        "length": length,
        "N_X": n_x,
        "N_rot": n_rot,
        "reduced_length": len(freely_reduce(word)),
    }
    if expected is not None:
        info["matches_expected"] = (result == tuple(expected))
        return info["matches_expected"], info
    return True, info


def sorts(word, state):
    """True iff the word sorts `state` to the identity."""
    n = len(state)
    ok, info = verify_word(tuple(state), word, tuple(range(n)))
    return ok


if __name__ == "__main__":
    from moves import identity, sigma
    # geodesic examples from PROBLEM.md §4.2 lead from id_n to sigma_n
    examples = {
        5: "XRXRRXRRXR",
        6: "XRXLXLLXLXRXRRR",
        7: "XLXRXRXLXLXLLXLXRXRRR",
        8: "XRXLXLXRXRXRRXRXRXLXLXRXRRRR",
    }
    for n, w in examples.items():
        ok, info = verify_word(identity(n), w, sigma(n))
        print(f"n={n}: word of length {info['length']} "
              f"(N_X={info['N_X']}, N_rot={info['N_rot']}) "
              f"reaches sigma_n: {ok}")
        assert ok and info["length"] == n * (n - 1) // 2
    print("all PROBLEM §4.2 example words verified")
