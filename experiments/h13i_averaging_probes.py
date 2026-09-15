"""h13i_averaging_probes-1.0

Probes three averaging-style candidate proofs of H13-I
(`I(pi) <= floor((n-1)^2/4)` for every permutation pi of Z_n, `I(pi)` the
minimum number of inversions over all n^2 double cuts (q, c); see
`docs/notes/h13_line_model.md`). All three are exhaustive refutations at
small n, kept together with the exact "sum over all n^2 origins" formula
used to state them precisely.

Definitions match `docs/notes/h13_line_model.md` §0: for a cut with position
start s and value origin v0, the line is
b_m = pi[(s+m) mod n], w_m = (b_m - v0) mod n, m = 0..n-1;
inv(s, v0) = number of inversions of w. I(pi) = min_{s,v0} inv(s, v0).

Probe A (joint average): min_{s,v0} inv(s, v0) <= (1/n^2) sum_{s,v0} inv(s, v0).
  `formula_avg_num` computes n^2 * (this average) in closed form (no loop
  over s, v0) via the identity, for a fixed unordered pair of positions
  {i, i'} with d = (i'-i) mod n, e = (pi(i')-pi(i)) mod n (both in [1, n-1]):
  the pair is inverted for exactly n(d+e) - 2de of the n^2 origins.
  `check_formula_matches_bruteforce` certifies the identity against a
  literal double loop over s, v0.

Probe B (exact-min over v0, then average over s): for each s, min_v0
inv(s, v0) is computed exactly (inner minimum, not an average); the
candidate bound is min_s [that exact minimum].

Probe C (datapoint pivots): restrict origins to the n cuts (s, v0) = (i0,
pi(i0)) for i0 = 0..n-1 (i.e. the cut is anchored at one of the pi points
itself), and take the exact minimum of inv over just these n origins.

All three are refuted (there is a pi with the probed quantity > floor((n-1)^2/4))
already at n = 4. Probes A/B/C are candidate *upper bounds* for I(pi), not
alternative definitions; I(pi) itself (the true minimum over all n^2 origins)
is separately VERIFIED <= floor((n-1)^2/4) for 4 <= n <= 10 (C33).
"""
import argparse
import itertools
import json
import time


def floor_sq(n):
    return ((n - 1) ** 2) // 4


def inv_of(pi, s, v0, n):
    b = [(pi[(s + m) % n] - v0) % n for m in range(n)]
    return sum(1 for j in range(n) for k in range(j + 1, n) if b[j] > b[k])


def exact_I(pi, n):
    best = None
    for s in range(n):
        for v0 in range(n):
            v = inv_of(pi, s, v0, n)
            if best is None or v < best:
                best = v
    return best


def formula_avg_num(pi, n):
    """n^2 * (average of inv over all n^2 origins), closed form."""
    total = 0
    for i in range(n):
        for ip in range(i + 1, n):
            d = (ip - i) % n
            e = (pi[ip] - pi[i]) % n
            total += n * (d + e) - 2 * d * e
    return total


def probe_A_num(pi, n):
    return formula_avg_num(pi, n)  # candidate bound: probe_A_num / n^2


def probe_B_num(pi, n):
    """min over s of (sum over v0 of inv(s, v0)); candidate bound: this / n."""
    best = None
    for s in range(n):
        total = 0
        b = [pi[(s + m) % n] for m in range(n)]
        for j in range(n):
            for k in range(j + 1, n):
                total += (b[k] - b[j]) % n
        if best is None or total < best:
            best = total
    return best


def probe_C(pi, n):
    """exact min of inv over the n datapoint-anchored origins (s, v0) = (i0, pi[i0])."""
    best = None
    for i0 in range(n):
        v = inv_of(pi, i0, pi[i0], n)
        if best is None or v < best:
            best = v
    return best


def check_formula_matches_bruteforce(nmax):
    for n in range(4, nmax + 1):
        sample = itertools.permutations(range(n))
        checked = 0
        for pi in sample:
            brute_total = sum(
                inv_of(pi, s, v0, n) for s in range(n) for v0 in range(n)
            )
            assert formula_avg_num(pi, n) == brute_total, (n, pi)
            checked += 1
            if checked >= 30:
                break
    return True


def run(nmax):
    check_formula_matches_bruteforce(min(nmax, 7))
    results = []
    for n in range(4, nmax + 1):
        bound = floor_sq(n)
        t0 = time.time()
        worst_A = worst_B = worst_C = -1
        worst_A_pi = worst_B_pi = worst_C_pi = None
        count = 0
        for pi in itertools.permutations(range(n)):
            a = probe_A_num(pi, n)  # compare a > bound * n^2
            if a > worst_A:
                worst_A, worst_A_pi = a, pi
            b = probe_B_num(pi, n)  # compare b > bound * n
            if b > worst_B:
                worst_B, worst_B_pi = b, pi
            c = probe_C(pi, n)  # compare c > bound
            if c > worst_C:
                worst_C, worst_C_pi = c, pi
            count += 1
        dt = time.time() - t0
        row = dict(
            n=n,
            bound=bound,
            perms_checked=count,
            probe_A_ok=worst_A <= bound * n * n,
            probe_A_worst_avg=worst_A / (n * n),
            probe_A_worst_pi=worst_A_pi,
            probe_B_ok=worst_B <= bound * n,
            probe_B_worst_avg=worst_B / n,
            probe_B_worst_pi=worst_B_pi,
            probe_C_ok=worst_C <= bound,
            probe_C_worst=worst_C,
            probe_C_worst_pi=worst_C_pi,
            seconds=round(dt, 2),
        )
        results.append(row)
        print(
            f"n={n} bound={bound} perms={count} "
            f"A: worst_avg={row['probe_A_worst_avg']:.3f} ok={row['probe_A_ok']} "
            f"B: worst_avg={row['probe_B_worst_avg']:.3f} ok={row['probe_B_ok']} "
            f"C: worst={worst_C} ok={row['probe_C_ok']} ({dt:.1f}s)"
        )
    return results


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--json-out", type=str, default=None)
    args = ap.parse_args()
    results = run(args.nmax)
    if args.json_out:
        with open(args.json_out, "w") as f:
            json.dump(results, f, indent=2)
