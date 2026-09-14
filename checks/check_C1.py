"""Independent checks of C1 (docs/incoming/lrx_2n3_proof_checked.md), the
baseline 2n^2/3-flavoured upper bound `D_n <= floor(A_n)`.

Written independently from the text after a logical audit (see
docs/notes/audit_C1.md): the cycle-load bookkeeping (loads h_i, M, t, E),
the carrier construction (Section 2: contraction, XL/RX macros, disjoint
XX-cancellation), the per-shift bound (Section 4: f_c, F_c, N_c, H(c)),
and the averaging identities (Section 5, eq. 11-12) are all reimplemented
from the proof text, not copied from docs/incoming/lrx_2n3_check.py (which
is kept, unmodified, as a second independent construction to cross-check
lengths against — AGENTS.md rule 3: an independent audit may re-implement
the moves from the specification). Every constructed word is applied
letter by letter with the reference moves of oracle/moves.py (oracle-1.0).

Checks, each tagged with the equation/section it targets:
- eq(1): the load step h_i - h_{i-1} at cycle vertices, by direct
  recount of edge loads (not by trusting the recurrence).
- eq(2): F + k <= n*M + E, for every nontrivial cycle.
- Section 1 carrier existence: a vertex b with h_b = M and a -to+ sign
  change exists whenever signs are mixed.
- eq(4): the constructed internal word W_C has length <= 2F - S(C),
  S(C) = 2M + E - 1, and actually undoes the cycle (index permutation
  C^{-1}) while fixing every slot outside C -- checked by application.
- eq(5): 3n*S(C) >= 5(F+k), for every nontrivial cycle (the case-split
  lemma of Section 3).
- eq(9): the universal route visits every vertex and has length H(c).
- eq(10): the full per-shift word sorts pi and has length within the
  stated per-shift bound.
- eq(11)-(12): the averaging identities, exact (Fraction), for sampled pi.
- The final theorem: constructed length <= floor(A_n) = baseline_bound_2n3(n)
  (bounds/known.py), exhaustively for small n and sampled for larger n;
  cross-checked against the certified true distance d(pi) where available
  (n <= 12) and against docs/incoming/lrx_2n3_check.py's own construction.
- Section 6: the sharpness family C_r (n=4r) achieves eq(5) with equality.

Passing finite checks do not prove the general theorem (AGENTS.md rule 5);
they test the implementation and hunt for errors in the lemmas on small
inputs, and, cross-checked against real BFS tables, show how loose the
bound is in practice.

Usage: python checks/check_C1.py [--cycle-max-n 8] [--perm-max-n 7]
       [--table-max-n 12] [--large 20,21,100,101] [--quick]
Output: data/runs/check_C1/report.json and summary.md
Version check_C1-1.0.
"""

import argparse
import itertools
import json
import os
import random
import sys
import time
from fractions import Fraction

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
sys.path.insert(0, os.path.join(ROOT, "exact"))
sys.path.insert(0, os.path.join(ROOT, "bounds"))
sys.path.insert(0, os.path.join(ROOT, "docs", "incoming"))

from moves import apply_word, freely_reduce, identity, delta, CORE_VERSION  # noqa: E402
from bfs import rank_perm, factorials  # noqa: E402
import known  # noqa: E402
import lrx_2n3_check as ref  # noqa: E402  (the input script -- independent second construction)

OUT_DIR = os.path.join(ROOT, "data", "runs", "check_C1")
TABLES_DIR = os.path.join(ROOT, "data", "tables")


class CheckFailure(Exception):
    pass


def check(cond, tag, info=None):
    if not cond:
        raise CheckFailure(f"{tag}: {info!r}")


# ---------------------------------------------------------------- Section 1

def signed_shortest_step(a, b, n):
    """Signed length of a shortest path a -> b on Z_n; positive direction
    chosen at an antipodal tie (n even, |b-a| = n/2) -- an explicit instance
    of the text's "choose either direction in a tie"."""
    fwd = (b - a) % n
    bwd = fwd - n
    return fwd if fwd <= -bwd else bwd


def path_edges(a, d, n):
    """Edge indices used by the signed path of length d starting at a;
    edge i means {i, i+1}."""
    if d > 0:
        return [(a + j) % n for j in range(d)]
    return [(a - j - 1) % n for j in range(-d)]


def cycle_stats(cycle, n):
    """F, per-edge loads, M, t, E, and the signed steps d_j for a cycle."""
    k = len(cycle)
    ds = [signed_shortest_step(cycle[j], cycle[(j + 1) % k], n) for j in range(k)]
    loads = [0] * n
    for a, d in zip(cycle, ds):
        for e in path_edges(a, d, n):
            loads[e] += 1
    F = sum(abs(d) for d in ds)
    M = max(loads)
    t = sum(1 for j in range(k) if ds[j - 1] < 0 < ds[j])
    E = k - 2 * t
    return ds, loads, F, M, t, E


def check_eq1_load_step(cycle, n, ds, loads):
    """eq(1): h_i - h_{i-1} at each cycle vertex i = a_j matches the sign
    pattern of (d_{j-1}, d_j); checked by direct recount, not by trusting
    the recurrence to build `loads`."""
    k = len(cycle)
    for j in range(k):
        i = cycle[j]
        h_i = loads[i % n]
        h_im1 = loads[(i - 1) % n]
        dprev, dcur = ds[j - 1], ds[j]
        if dprev < 0 < dcur:
            expect = 2
        elif dprev > 0 > dcur:
            expect = -2
        else:
            expect = 0
        check(h_i - h_im1 == expect, "eq1", (cycle, i, dprev, dcur, h_i, h_im1))


def find_carrier(cycle, ds, loads, M, n):
    """Section 1: a vertex b = a_j with h_b = M and a -to+ sign change,
    when signs are mixed; any vertex, when all signs agree."""
    k = len(cycle)
    mixed = any(d > 0 for d in ds) and any(d < 0 for d in ds)
    if not mixed:
        return 0
    for j in range(k):
        if ds[j - 1] < 0 < ds[j] and loads[cycle[j] % n] == M:
            return j
    raise CheckFailure(f"no -to+ max-load carrier found: cycle={cycle}")


# ---------------------------------------------------------------- Section 2

def build_cycle_word(cycle, n):
    """The internal word W_C (Section 2): rotate the cycle to start at the
    carrier, treat it as fixed, and walk the remaining ring with XL/RX
    macros (contracting the edge at the carrier), doing one X per exchange.
    Returns (word, carrier_position, F, M, t, E)."""
    k = len(cycle)
    ds, loads, F, M, t, E = cycle_stats(cycle, n)
    start = find_carrier(cycle, ds, loads, M, n)
    check_eq1_load_step(cycle, n, ds, loads)

    rot = cycle[start:] + cycle[:start]
    carrier = rot[0]
    rds = [signed_shortest_step(rot[j], rot[(j + 1) % k], n) for j in range(k)]
    contracted = carrier

    raw = []
    for j in range(k):
        a, d = rot[j], rds[j]
        steps = sum(1 for e in path_edges(a, d, n) if e != contracted)
        macro = ("X", "L") if d > 0 else ("R", "X")
        raw.extend(macro * steps)
        if j < k - 1:
            raw.append("X")
    word = freely_reduce("".join(raw))

    S = 2 * M + E - 1
    check(len(word) <= 2 * F - S, "eq4-length", (cycle, len(word), F, S))
    return word, carrier, F, M, t, E


def check_eq2(F, k, n, M, E):
    check(F + k <= n * M + E, "eq2", (F, k, n, M, E))


def check_eq5(F, k, n, M, E):
    S = 2 * M + E - 1
    check(3 * n * S >= 5 * (F + k), "eq5", (F, k, n, M, E, S))


def check_cycle_word_action(cycle, carrier, word, n):
    """W_C is defined (Section 2) relative to a frame already rotated so
    that `carrier` sits at array index 0 ("the initial positioning
    rotation", explicitly excluded from W_C itself). So the isolated check
    is: L^carrier, then W_C, then R^carrier (undo the framing) realizes
    C^{-1} as the array-index permutation restricted to the cycle's
    slots, restoring the carrier and everything outside the cycle."""
    k = len(cycle)
    g = list(range(n))
    for j in range(k):
        g[cycle[j]] = cycle[j - 1]  # C^{-1}: successor -> predecessor
    target = tuple(g)
    framed = "L" * carrier + word + "R" * carrier
    got = apply_word(identity(n), framed)
    check(got == target, "cycle-word-action", (cycle, word, got, target))


# ---------------------------------------------------------------- Section 4

def universal_route(n, c):
    """H(c) route (eq 9): sequence of (vertex, move) visiting every vertex,
    ending at c."""
    c %= n
    route = []
    cur = 0
    if c == 0:
        for _ in range(n):
            cur = (cur + 1) % n
            route.append((cur, "L"))
        return route
    d = signed_shortest_step(0, c, n)
    short_move = "L" if d > 0 else "R"
    long_move = "R" if d > 0 else "L"
    for _ in range(n - 1):
        cur = (cur + (1 if long_move == "L" else -1)) % n
        route.append((cur, long_move))
    for _ in range(abs(d) - 1):
        cur = (cur + (1 if short_move == "L" else -1)) % n
        route.append((cur, short_move))
    check(cur == c, "eq9-endpoint", (n, c, cur))
    return route


def H_of_c(n, c):
    return n if c % n == 0 else n + delta(0, c % n, n) - 2


def nontrivial_cycles_of(mapping):
    n = len(mapping)
    seen = [False] * n
    cycles = []
    for s in range(n):
        if seen[s]:
            continue
        cyc = []
        x = s
        while not seen[x]:
            seen[x] = True
            cyc.append(x)
            x = mapping[x]
        if len(cyc) > 1:
            cycles.append(tuple(cyc))
    return cycles


def word_for_shift(state, c):
    n = len(state)
    f_c = [(v + c) % n for v in state]
    cycles = nontrivial_cycles_of(f_c)
    plans = {}
    F_c = sum(delta(i, f_c[i], n) for i in range(n))
    N_c = sum(1 for i in range(n) if f_c[i] != i)
    total_F = total_k = 0
    for cyc in cycles:
        word, carrier, F, M, t, E = build_cycle_word(list(cyc), n)
        check_eq2(F, len(cyc), n, M, E)
        check_eq5(F, len(cyc), n, M, E)
        check_cycle_word_action(list(cyc), carrier, word, n)
        plans[carrier] = word
        total_F += F
        total_k += len(cyc)
    check(total_F == F_c, "F_c-matches-cycle-sum", (state, c, total_F, F_c))

    raw = []
    used = set()
    if 0 in plans:
        raw.append(plans[0])
        used.add(0)
    route = universal_route(n, c)
    check(len(route) <= H_of_c(n, c), "eq9-length", (n, c, len(route), H_of_c(n, c)))
    for vertex, move in route:
        raw.append(move)
        if vertex in plans and vertex not in used:
            raw.append(plans[vertex])
            used.add(vertex)
    check(used == set(plans), "route-visits-all-carriers", (state, c, used, set(plans)))
    word = freely_reduce("".join(raw))
    return word, F_c, N_c


def check_eq10(state, c, word, F_c, N_c, n):
    """d(pi,id) <= (2 - 5/(3n))*F_c - 5/(3n)*N_c + H(c) (eq 10); checked
    against the ACTUAL constructed word length (a stronger, sufficient
    fact -- the cycles' own eq(4)+eq(5) bounds sum to this)."""
    bound = Fraction(2, 1) * F_c - Fraction(5, 3 * n) * F_c - Fraction(5, 3 * n) * N_c + H_of_c(n, c)
    check(len(word) <= bound, "eq10", (state, c, len(word), bound))


# ---------------------------------------------------------------- Section 5

def check_averaging_identities(n, samples):
    """eq(11)-(12): identities that hold for EVERY pi (F_c, N_c) or are
    pi-independent (H(c)); exact via Fraction."""
    P = n * n // 4
    for state in samples:
        sum_F = sum(sum(delta(i, (state[i] + c) % n, n) for i in range(n)) for c in range(n))
        sum_N = sum(sum(1 for i in range(n) if (state[i] + c) % n != i) for c in range(n))
        check(Fraction(sum_F, n) == P, "eq11-F", (state, sum_F, n, P))
        check(Fraction(sum_N, n) == n - 1, "eq11-N", (state, sum_N, n))
    sum_H = sum(H_of_c(n, c) for c in range(n))
    check(Fraction(sum_H, n) == n - 2 + Fraction(P + 2, n), "eq12-H", (n, sum_H))


# ---------------------------------------------------------------- top level

def sort_and_bound(state):
    n = len(state)
    best = None
    for c in range(n):
        word, F_c, N_c = word_for_shift(state, c)
        check_eq10(state, c, word, F_c, N_c, n)
        got = apply_word(state, word)
        check(got == identity(n), "sorts-to-identity", (state, c, word, got))
        if best is None or (len(word), c) < (len(best[0]), best[1]):
            best = (word, c)
    return best


def cycles_of_zn(n, max_full=None):
    """All nontrivial cycles (as position-sequences) of length 2..n on Z_n,
    i.e. all injective sequences up to rotation (fixing the first element).
    max_full caps enumeration for large n (sampled instead)."""
    out = []
    for k in range(2, n + 1):
        rest = [x for x in range(n) if x != 0]
        cnt = 0
        for perm in itertools.permutations(rest, k - 1):
            out.append((0,) + perm)
            cnt += 1
            if max_full is not None and cnt >= max_full:
                break
    return out


def run(cycle_max_n, perm_max_n, table_max_n, large_ns, quick, seed=1):
    rng = random.Random(seed)
    t0 = time.time()
    report = {
        "version": "check_C1-1.0",
        "core_version": CORE_VERSION,
        "source": "docs/incoming/lrx_2n3_proof_checked.md",
        "cycles": {},
        "permutations": {},
        "tables": {},
        "large": {},
        "sharpness_family": {},
        "averaging": {},
        "scalar": {},
    }

    # --- cycles: eq(1),(2),(4),(5), action ---
    # cycles_of_zn caps enumeration PER cycle length k at max_full; the
    # largest single-k count for n <= 9 is (n-1)! at k=n, which stays under
    # this cap, so every n in the default range (<=8) is exhaustive over
    # ALL nontrivial cycles of Z_n (all lengths, all starting orientations).
    cap = 50000
    for n in range(4, cycle_max_n + 1):
        largest_k_count = 1
        for j in range(1, n):  # (n-1)!, the count at k=n (full cycle)
            largest_k_count *= j
        exhaustive = largest_k_count <= cap
        cycles = cycles_of_zn(n, max_full=cap)
        checked = 0
        for cyc in cycles:
            word, carrier, F, M, t, E = build_cycle_word(list(cyc), n)
            check_eq2(F, len(cyc), n, M, E)
            check_eq5(F, len(cyc), n, M, E)
            check_cycle_word_action(list(cyc), carrier, word, n)
            checked += 1
        report["cycles"][n] = {"checked": checked, "coverage": "exhaustive" if exhaustive else f"sampled(cap={cap} per cycle length)"}
        print(f"cycles n={n}: {checked} checked, PASS")

    # --- averaging identities ---
    for n in range(4, 9):
        samples = [tuple(range(n))]
        for _ in range(5):
            p = list(range(n))
            rng.shuffle(p)
            samples.append(tuple(p))
        check_averaging_identities(n, samples)
        report["averaging"][n] = {"samples": len(samples)}
    print("averaging identities eq(11)-(12): PASS")

    # --- full permutations, small n: exhaustive sort + bound ---
    for n in range(4, perm_max_n + 1):
        bound = known.baseline_bound_2n3(n)
        maxlen = -1
        checked = 0
        for state in itertools.permutations(range(n)):
            word, c = sort_and_bound(state)
            check(len(word) <= bound, "theorem-bound", (state, len(word), bound))
            maxlen = max(maxlen, len(word))
            checked += 1
        report["permutations"][n] = {"checked": checked, "bound": bound, "max_constructed": maxlen}
        print(f"permutations n={n}: {checked} checked, bound={bound}, max_constructed={maxlen}, PASS")

    # --- cross-check vs certified true distances, n up to table_max_n ---
    fact_cache = {}
    for n in list(range(perm_max_n + 1, table_max_n + 1)):
        path = os.path.join(TABLES_DIR, f"dist_n{n}.bin")
        if not os.path.exists(path):
            continue
        fact = factorials(n)
        with open(path, "rb") as f:
            dist = f.read()
        bound = known.baseline_bound_2n3(n)
        n_samples = 30 if not quick else 5
        worst_gap = -1
        for _ in range(n_samples):
            p = list(range(n))
            rng.shuffle(p)
            state = tuple(p)
            word, c = sort_and_bound(state)
            check(len(word) <= bound, "theorem-bound", (state, len(word), bound))
            true_d = dist[rank_perm(state, fact)]
            check(len(word) >= true_d, "constructed-not-below-true-d", (state, len(word), true_d))
            worst_gap = max(worst_gap, len(word) - true_d)
        report["tables"][n] = {"samples": n_samples, "bound": bound, "max_gap_over_true_d": worst_gap}
        print(f"n={n} vs certified table: {n_samples} samples, bound={bound}, max(len-d)={worst_gap}, PASS")

    # --- large n: bound only (+ ref-script cross-check) ---
    for n in large_ns:
        bound = known.baseline_bound_2n3(n)
        n_samples = 10 if not quick else 2
        maxlen = -1
        maxlen_ref = -1
        for _ in range(n_samples):
            p = list(range(n))
            rng.shuffle(p)
            state = tuple(p)
            word, c = sort_and_bound(state)
            check(len(word) <= bound, "theorem-bound", (state, len(word), bound))
            ref_result = ref.sorting_word(state)
            ref_final = ref.apply_word(state, ref_result.word)
            check(ref_final == tuple(range(n)), "ref-script-sorts", (state, ref_result.word))
            check(len(ref_result.word) <= bound, "ref-script-bound", (state, len(ref_result.word), bound))
            maxlen = max(maxlen, len(word))
            maxlen_ref = max(maxlen_ref, len(ref_result.word))
        report["large"][n] = {
            "samples": n_samples, "bound": bound,
            "max_constructed": maxlen, "max_constructed_ref_script": maxlen_ref,
        }
        print(f"n={n}: {n_samples} samples, bound={bound}, "
              f"max_len(mine)={maxlen}, max_len(ref)={maxlen_ref}, PASS")

    # --- Section 6: sharpness family C_r, n=4r ---
    for r in range(2, 9):
        n = 4 * r
        cyc = [0, 2 * r - 1, r, 3 * r - 1, 2 * r, 4 * r - 1, 3 * r, r - 1]
        ds, loads, F, M, t, E = cycle_stats(cyc, n)
        k = 8
        check(t == 4 and E == 0 and F == 3 * n - 8 and M == 3, "sharpness-stats",
              (r, n, t, E, F, M))
        S = 2 * M + E - 1
        check(3 * n * S == 5 * (F + k), "sharpness-equality", (r, n, S, F, k))
        report["sharpness_family"][r] = {"n": n, "F": F, "M": M, "t": t, "E": E, "S": S}
    print("sharpness family C_r (Section 6): PASS")

    # --- scalar: independent ceiling formula vs bounds/known.py, and 2nd bound ---
    mismatches = []
    for n in list(range(4, 2001)) + [5000, 10000, 100000]:
        P = n * n // 4
        num = 2 * P + 11 * n - 11
        ceiling = -(-num // (3 * n))
        mine = 2 * P + n - ceiling
        theirs = known.baseline_bound_2n3(n)
        if mine != theirs:
            mismatches.append((n, mine, theirs))
        second_bound = n * (n - 1) // 2 + (4 * n // 3) - 3
        if n >= 5:
            check(theirs <= second_bound, "second-bound-dominates", (n, theirs, second_bound))
        elif n == 4:
            check(theirs == second_bound == 8, "n4-bounds-equal", (n, theirs, second_bound))
    check(not mismatches, "scalar-formula-matches-known", mismatches[:5])
    report["scalar"] = {"n_checked": 2000 + 3, "mismatches": len(mismatches)}
    print(f"scalar formula cross-check (independent ceiling vs bounds/known.py): PASS, "
          f"{report['scalar']['n_checked']} values of n")

    report["elapsed_s"] = time.time() - t0
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(os.path.join(OUT_DIR, "report.json"), "w") as f:
        json.dump(report, f, indent=2, default=str)
    with open(os.path.join(OUT_DIR, "summary.md"), "w") as f:
        f.write(f"# check_C1 summary\n\nversion {report['version']}, core {CORE_VERSION}\n\n")
        f.write(f"elapsed: {report['elapsed_s']:.1f} s\n\n")
        f.write("All checks PASS. See report.json for coverage detail.\n")
    print(f"\nALL CHECKS PASS ({report['elapsed_s']:.1f} s). Report: {OUT_DIR}")
    return report


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cycle-max-n", type=int, default=8)
    ap.add_argument("--perm-max-n", type=int, default=7)
    ap.add_argument("--table-max-n", type=int, default=12)
    ap.add_argument("--large", type=str, default="20,21,100,101")
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    large_ns = [int(x) for x in args.large.split(",") if x]
    run(args.cycle_max_n, args.perm_max_n, args.table_max_n, large_ns, args.quick)
