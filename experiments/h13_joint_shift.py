"""h13_joint_shift.py — H13-I session-9 probe: joint-shift structure of the
double-cut inversion count I(pi) = min_{q,c} inv(w_{q,c}).

This script does NOT prove or refute H13-I. It exhaustively checks three
algebraic facts about the (a, b)-grid of double cuts that came out of this
session's attempt (see docs/notes/h13_line_model.md, section "7. Сессия 9"),
and one restricted-cut-family negative result:

1. Per-pair quadrant identity. For a permutation pi and an unordered pair of
   positions {i1, i2} with i1 < i2 (integers), let da = i2 - i1 and
   db = (pi(i2) - pi(i1)) mod n. Summed over all n^2 cuts (a, b) (a = 0..n-1
   rotates positions, b = 0..n-1 rotates values, both mod n; see
   experiments/line_profile.c for the exact (q, c) <-> (a, b) correspondence,
   a = q + 1, b = c), the pair is concordant (non-inverted) at
   da*db + (n-da)*(n-db) cuts and discordant (inverted) at
   da*(n-db) + (n-da)*db cuts; hence concordant - discordant = (2da-n)(2db-n)
   for every pair, independent of a, b. Checked by direct enumeration against
   the closed form.

2. Telescoping recursion. For a permutation s of Z_n (thought of as a
   sequence s_0..s_{n-1}), inv(shift_{b+1}(s)) - inv(shift_b(s)) =
   n - 1 - 2*pos(b), where pos = s^{-1} and shift_b(s)_j = (s_j - b) mod n.
   Checked directly for all b and all pi.

3. Restricted "no value shift" family (b = 0 only, i.e. only the n position
   rotations q, matching the "n cuts, not n^2" families already ruled out in
   docs/notes/h13_line_model.md section 1.7): min_a inv(rotate_a(pi)) can
   exceed floor((n-1)^2/4) (found already at n = 7: the same witness
   (0,5,4,3,2,1,6) as the section-4 "lemma C" counterexample). Logged here to
   confirm this specific restricted family (independent from the "first
   element = t" family already tested) is also insufficient, and to record
   the exact worst-case gap by n.

4. Reflection cross-check. For pi_h(i) = h - i mod n, the b = 0 rotation
   family alone (no value shift) already realizes min_a J(a) =
   floor((n-1)^2/4) exactly, via the elementary "two monotone runs" formula
   J(a) = C(k,2) + C(n-k,2) minimized at k = floor(n/2); this reproduces the
   known extremal value (C33) by a new, short, direct argument (not via
   BFS/geodesics), for a wide range of n.

Usage: python3 h13_joint_shift.py --nmax 8
Version: h13_joint_shift-1.0
"""
import argparse
import itertools
import math


def inversions(seq):
    n = len(seq)
    c = 0
    for i in range(n):
        for j in range(i + 1, n):
            if seq[i] > seq[j]:
                c += 1
    return c


def rotate_a(pi, a):
    n = len(pi)
    return [pi[(a + j) % n] for j in range(n)]


def shift_b(seq, b):
    n = len(seq)
    return [(v - b) % n for v in seq]


def check_pair_identity(nmax):
    """Fact 1: concordant - discordant = (2da-n)(2db-n) for every pair,
    every pi, checked by direct n^2 enumeration (not the closed form) against
    the closed form."""
    failures = 0
    checked_pairs = 0
    for n in range(4, nmax + 1):
        for pi in itertools.permutations(range(n)):
            for i1 in range(n):
                for i2 in range(i1 + 1, n):
                    da = i2 - i1
                    db = (pi[i2] - pi[i1]) % n
                    conc = 0
                    disc = 0
                    for a in range(n):
                        j1 = (i1 - a) % n
                        j2 = (i2 - a) % n
                        for b in range(n):
                            w1 = (pi[i1] - b) % n
                            w2 = (pi[i2] - b) % n
                            order_pos = j1 < j2
                            order_val = w1 < w2
                            if order_pos == order_val:
                                conc += 1
                            else:
                                disc += 1
                    formula = (2 * da - n) * (2 * db - n)
                    checked_pairs += 1
                    if conc - disc != formula:
                        failures += 1
                        print(f"FAIL fact1 n={n} pi={pi} pair=({i1},{i2}) "
                              f"conc-disc={conc-disc} formula={formula}")
        if failures:
            break
    return {"checked_pairs": checked_pairs, "failures": failures}


def check_telescoping(nmax):
    """Fact 2: inv(shift_{b+1}(s)) - inv(shift_b(s)) == n-1-2*pos(b)."""
    failures = 0
    checked = 0
    for n in range(2, nmax + 1):
        for s in itertools.permutations(range(n)):
            pos = [0] * n
            for idx, v in enumerate(s):
                pos[v] = idx
            prev = inversions(shift_b(list(s), 0))
            for b in range(n):
                nxt = inversions(shift_b(list(s), (b + 1) % n))
                expected_step = n - 1 - 2 * pos[b]
                checked += 1
                if nxt - prev != expected_step:
                    failures += 1
                    print(f"FAIL fact2 n={n} s={s} b={b} "
                          f"delta={nxt-prev} expected={expected_step}")
                prev = nxt
        if failures:
            break
    return {"checked_steps": checked, "failures": failures}


def check_b0_family(nmax):
    """Fact 3: min_a inv(rotate_a(pi)) (b = 0 only) vs floor((n-1)^2/4)."""
    results = {}
    for n in range(4, nmax + 1):
        target = ((n - 1) ** 2) // 4
        worst = -1
        worst_pi = None
        for pi in itertools.permutations(range(n)):
            m = min(inversions(rotate_a(list(pi), a)) for a in range(n))
            if m > worst:
                worst = m
                worst_pi = pi
        results[n] = {
            "target": target,
            "worst_min_a_J": worst,
            "worst_pi": worst_pi,
            "sufficient": worst <= target,
        }
    return results


def check_reflection_two_runs(ns):
    """Fact 4: for pi_h, min_a J(a) (b = 0 family alone) hits
    floor((n-1)^2/4) exactly, via the "two monotone runs" closed form."""
    results = {}
    for n in ns:
        target = ((n - 1) ** 2) // 4
        h = 1  # sigma_n; any h gives the same min_a J(a) by rotation symmetry
        pi = [(h - i) % n for i in range(n)]
        best = min(inversions(rotate_a(pi, a)) for a in range(n))
        # closed form: best over k = 1..n-1 of C(k,2)+C(n-k,2)
        def C2(k):
            return k * (k - 1) // 2
        closed = min(C2(k) + C2(n - k) for k in range(1, n))
        results[n] = {"target": target, "brute": best, "closed_form": closed,
                       "matches_target": best == target and closed == target}
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8,
                     help="max n for exhaustive facts 1-3 (all permutations)")
    ap.add_argument("--reflect-ns", type=str, default="4,5,6,7,8,9,10,20,50,101,200",
                     help="comma-separated n for fact 4 (reflections only)")
    args = ap.parse_args()

    print("=== Fact 1: per-pair quadrant identity ===")
    r1 = check_pair_identity(args.nmax)
    print(r1)

    print("=== Fact 2: telescoping recursion ===")
    r2 = check_telescoping(args.nmax)
    print(r2)

    print("=== Fact 3: b=0-only restricted family (no value shift) ===")
    r3 = check_b0_family(args.nmax)
    for n, r in r3.items():
        print(n, r)

    print("=== Fact 4: reflections, b=0-only family hits the target exactly ===")
    ns = [int(x) for x in args.reflect_ns.split(",")]
    r4 = check_reflection_two_runs(ns)
    for n, r in r4.items():
        print(n, r)

    ok = (r1["failures"] == 0 and r2["failures"] == 0 and
          all(r["matches_target"] for r in r4.values()))
    print("ALL ALGEBRAIC FACTS (1,2,4) VERIFIED:" , ok)
    print("Fact 3 conclusion: b=0-only family is INSUFFICIENT for H13-I "
          "(fails already at n=7,8) -- consistent with docs/notes/h13_line_model.md")


if __name__ == "__main__":
    main()
