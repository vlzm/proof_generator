#!/usr/bin/env python3
"""Construct and verify the L/R/X sorting words from the edge-load proof.

The state is a tuple ``p`` in one-line notation: position i contains p[i].
Moves act from left to right:

    L [x0,x1,...,x(n-1)] = [x1,...,x(n-1),x0]
    R = L^{-1}
    X [x0,x1,...]       = [x1,x0,...]

For every c in Z/nZ the program constructs the proof's word for
f_c(i) = p[i] + c (mod n), and returns the shortest of these n words.
"""

from __future__ import annotations

import argparse
import itertools
from dataclasses import dataclass
from typing import Iterable, Sequence


Move = str


def apply_move(state: tuple[int, ...], move: Move) -> tuple[int, ...]:
    """Apply one of L, R, X to a state."""
    if move == "L":
        return state[1:] + state[:1]
    if move == "R":
        return state[-1:] + state[:-1]
    if move == "X":
        return (state[1], state[0]) + state[2:]
    raise ValueError(f"unknown move {move!r}")


def apply_word(state: Sequence[int], word: Iterable[Move]) -> tuple[int, ...]:
    """Apply a word from left to right."""
    ans = tuple(state)
    for move in word:
        ans = apply_move(ans, move)
    return ans


def freely_reduce(word: Iterable[Move]) -> list[Move]:
    """Cancel XX, LR and RL without changing the represented permutation."""
    inverse = {"L": "R", "R": "L", "X": "X"}
    stack: list[Move] = []
    for move in word:
        if stack and stack[-1] == inverse[move]:
            stack.pop()
        else:
            stack.append(move)
    return stack


def signed_shortest_step(a: int, b: int, n: int) -> int:
    """Signed length of a shortest path a -> b on Z/nZ.

    At an even-n antipodal tie, the positive direction is chosen.
    """
    positive = (b - a) % n
    negative = positive - n
    return positive if positive <= -negative else negative


def path_edges(a: int, d: int, n: int) -> list[int]:
    """Edges used by a signed path, with edge i meaning {i,i+1}."""
    if d > 0:
        return [(a + j) % n for j in range(d)]
    return [(a - j - 1) % n for j in range(-d)]


def nontrivial_cycles(mapping: Sequence[int]) -> list[list[int]]:
    """Return the non-singleton cycles of a permutation mapping."""
    n = len(mapping)
    seen = [False] * n
    cycles: list[list[int]] = []
    for start in range(n):
        if seen[start]:
            continue
        cycle: list[int] = []
        x = start
        while not seen[x]:
            seen[x] = True
            cycle.append(x)
            x = mapping[x]
        if len(cycle) > 1:
            cycles.append(cycle)
    return cycles


@dataclass(frozen=True)
class CyclePlan:
    carrier: int
    word: tuple[Move, ...]
    cycle: tuple[int, ...]
    F: int
    M: int
    t: int


def plan_cycle(cycle: Sequence[int], n: int) -> CyclePlan:
    """Construct the internal word P_C from Sections 1--3 of the proof."""
    k = len(cycle)
    if k < 2:
        raise ValueError("plan_cycle expects a nontrivial cycle")

    ds = [
        signed_shortest_step(cycle[j], cycle[(j + 1) % k], n)
        for j in range(k)
    ]
    loads = [0] * n
    for a, d in zip(cycle, ds):
        for edge in path_edges(a, d, n):
            loads[edge] += 1

    M = max(loads)
    mixed = any(d > 0 for d in ds) and any(d < 0 for d in ds)
    if mixed:
        candidates = [
            j
            for j in range(k)
            if ds[j - 1] < 0 < ds[j] and loads[cycle[j]] == M
        ]
        if not candidates:
            raise AssertionError(
                "the maximum-load (- to +) carrier promised by the proof "
                "was not found"
            )
        start = candidates[0]
    else:
        # In the one-sign case the load is constant, so any vertex works.
        start = 0

    rotated = list(cycle[start:]) + list(cycle[:start])
    carrier = rotated[0]
    rds = [
        signed_shortest_step(rotated[j], rotated[(j + 1) % k], n)
        for j in range(k)
    ]

    raw: list[Move] = []
    contracted_edge = carrier  # edge (carrier, carrier+1)
    for j, (a, d) in enumerate(zip(rotated, rds)):
        steps = sum(edge != contracted_edge for edge in path_edges(a, d, n))
        # XL rotates the (n-1)-ring positively; RX rotates it negatively.
        macro = ("X", "L") if d > 0 else ("R", "X")
        raw.extend(macro * steps)
        if j < k - 1:
            raw.append("X")

    word = tuple(freely_reduce(raw))

    # This is both a programming invariant and a direct check of (3).
    t = sum(ds[j - 1] < 0 < ds[j] for j in range(k))
    E = k - 2 * t
    F = sum(abs(d) for d in ds)
    theoretical_cycle_limit = 2 * F - (2 * M + E - 1)
    if len(word) > theoretical_cycle_limit:
        raise AssertionError(
            f"cycle word has length {len(word)}, above (3)'s limit "
            f"{theoretical_cycle_limit}"
        )

    return CyclePlan(carrier, word, tuple(rotated), F, M, t)


def universal_route(n: int, c: int) -> list[tuple[int, Move]]:
    """The route in (12), returned as (new_vertex, move) unit steps."""
    c %= n
    route: list[tuple[int, Move]] = []
    current = 0
    if c == 0:
        for _ in range(n):
            current = (current + 1) % n
            route.append((current, "L"))
        return route

    d = signed_shortest_step(0, c, n)
    short_move = "L" if d > 0 else "R"
    long_move = "R" if d > 0 else "L"

    # Reach the neighbour of zero on c's side by traversing the long arc.
    for _ in range(n - 1):
        current = (current + (1 if long_move == "L" else -1)) % n
        route.append((current, long_move))
    # Retrace the short arc from that neighbour to c.
    for _ in range(abs(d) - 1):
        current = (current + (1 if short_move == "L" else -1)) % n
        route.append((current, short_move))
    if current != c:
        raise AssertionError("route construction did not end at c")
    return route


@dataclass(frozen=True)
class SortResult:
    word: tuple[Move, ...]
    shift: int


def word_for_shift(state: Sequence[int], c: int) -> tuple[Move, ...]:
    """Construct the complete word associated with a fixed final shift c."""
    n = len(state)
    f_c = [((value + c) % n) for value in state]
    plans = [plan_cycle(cycle, n) for cycle in nontrivial_cycles(f_c)]
    by_carrier = {plan.carrier: plan for plan in plans}
    if len(by_carrier) != len(plans):
        raise AssertionError("distinct cycles unexpectedly chose one carrier")

    raw: list[Move] = []
    used: set[int] = set()

    # A carrier at 0 is already available before the first route step.
    if 0 in by_carrier:
        raw.extend(by_carrier[0].word)
        used.add(0)

    for vertex, move in universal_route(n, c):
        raw.append(move)
        if vertex in by_carrier and vertex not in used:
            raw.extend(by_carrier[vertex].word)
            used.add(vertex)

    if len(used) != len(plans):
        raise AssertionError("universal route missed a cycle carrier")
    return tuple(freely_reduce(raw))


def sorting_word(state: Sequence[int]) -> SortResult:
    """Return the shortest proof-constructed word among all shifts c."""
    state = tuple(state)
    n = len(state)
    if sorted(state) != list(range(n)):
        raise ValueError("state must be a permutation of range(n)")
    candidates = [(word_for_shift(state, c), c) for c in range(n)]
    word, c = min(candidates, key=lambda item: (len(item[0]), item[1]))
    return SortResult(word, c)


def stated_bound(n: int) -> int:
    """The first (stronger) integer bound in equation (1)."""
    P = (n * n) // 4
    numerator = 2 * P + 11 * n - 11
    ceiling = (numerator + 3 * n - 1) // (3 * n)
    return 2 * P + n - ceiling


def exhaustive_check(n: int) -> dict[str, object]:
    """Execute every c-word, and the selected shortest word, on all of S_n."""
    identity = tuple(range(n))
    bound = stated_bound(n)
    maximum = -1
    worst_states: list[tuple[int, ...]] = []
    shift_counts = [0] * n
    checked = 0

    for state in itertools.permutations(range(n)):
        candidates: list[tuple[tuple[Move, ...], int]] = []
        for c in range(n):
            word_c = word_for_shift(state, c)
            final_c = apply_word(state, word_c)
            if final_c != identity:
                raise AssertionError(
                    f"n={n}: c-word failed for {state}; got {final_c}; "
                    f"c={c}, word={''.join(word_c)}"
                )
            candidates.append((word_c, c))

        word, shift = min(candidates, key=lambda item: (len(item[0]), item[1]))
        result = SortResult(word, shift)
        length = len(result.word)
        if length > bound:
            raise AssertionError(
                f"n={n}: bound {bound} failed for {state}; got {length}"
            )
        shift_counts[result.shift] += 1
        checked += 1
        if length > maximum:
            maximum = length
            worst_states = [state]
        elif length == maximum:
            worst_states.append(state)

    return {
        "n": n,
        "states": checked,
        "bound": bound,
        "maximum_constructed_length": maximum,
        "number_of_maximizers": len(worst_states),
        "first_maximizer": worst_states[0],
        "chosen_shift_counts": tuple(shift_counts),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--exhaustive",
        type=int,
        nargs="+",
        metavar="N",
        help="exhaustively verify all permutations for each supplied N",
    )
    parser.add_argument(
        "--state",
        type=int,
        nargs="+",
        help="sort one zero-based permutation, e.g. --state 2 0 1 3",
    )
    args = parser.parse_args()

    if not args.exhaustive and args.state is None:
        args.exhaustive = [4, 5, 6]

    if args.state is not None:
        result = sorting_word(args.state)
        final = apply_word(args.state, result.word)
        print(f"state:  {tuple(args.state)}")
        print(f"shift:  {result.shift}")
        print(f"length: {len(result.word)}")
        print(f"word:   {''.join(result.word) or '(empty)'}")
        print(f"final:  {final}")

    if args.exhaustive:
        for n in args.exhaustive:
            report = exhaustive_check(n)
            print(
                "n={n}: PASS; states={states}; bound={bound}; "
                "max constructed length={maximum_constructed_length}; "
                "maximizers={number_of_maximizers}; "
                "first maximizer={first_maximizer}; "
                "chosen c counts={chosen_shift_counts}".format(**report)
            )


if __name__ == "__main__":
    main()
