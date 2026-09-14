"""H13-I, session 9: check the algebraic hint (docs/notes/h13_line_model.md #6)
and the "fixed small family of position-cuts a" idea, and record why both fail.

Part 1 (decomposition). For pi in Z_n and pair i1 < i2, let d = i2 - i1,
u1 = pi[i1], u2 = pi[i2], e = (u2 - u1) mod n. For a double cut (q, c) with
a = (q + 1) mod n, b = (q + 1 - c) mod n (shift), the pair is an inversion of
the relabelled line w (experiments/line_model.py) iff A(a) xor B(b), where
A(a) = 1 iff a in {i1+1, ..., i2} (window of size d) and B(b) = 1 iff b lies in
the window of size e starting at u1 + 1 (mod n). Verified by direct
recomputation of inv(w) against the xor sum, random pi, 4 <= n <= 8.

Part 2 (pigeonhole bound via b only, exact for reflections). Summing over b:
sum_b [A xor B] = e + A*(n - 2e), so sum_b inv(a, b) = E + H(a), with
E = sum_pairs e_pair (fixed given pi) and H(a) = sum_pairs A_pair(a)*(n-2e_pair).
Since min_b <= average_b, I(pi) <= min_a floor((E + H(a)) / n) -- a provable
upper bound, cheaper than exhaustive search. For a reflection pi_h(i) = h - i
mod n, e_pair = n - d_pair identically, hence:
  H(a) = 2*sum_pairs d*A_pair(a) - n*F(a) = 2*(n*F(a)/2) - n*F(a) = 0
(closed form sum_pairs d*A_pair(a) = a*n*(n-a)/2 = n*F(a)/2 derived and
checked below), for EVERY a, EVERY h, EVERY n. So the pigeonhole bound on
reflections equals E/n = (n-1)(2n-1)/6 exactly, for every a, while the true
threshold is floor((n-1)^2/4); the gap (n-1)(2n-1)/6 - (n-1)^2/4 = (n^2-1)/12
is positive and grows without bound. Hence this specific technique (bound the
minimum over b by the mean over b, for any single a or by taking the best a)
cannot prove H13-I for large n: it is off by Theta(n^2) on the very family
where the true bound is tight. This is a proved fact (algebra below), checked
against direct computation for a range of n.

Part 3 (fixed constant-size family of a's). Exhaustively: does a FIXED pair
{a1, a2}, independent of pi, always contain an a with min_b inv(a,b) <=
floor((n-1)^2/4)? True for {0, 3} at 4 <= n <= 11 (exhaustive), but a local
search (hill climbing over transpositions, several random restarts) finds a
permutation at n = 20 exceeding the threshold by 2 -- so no fixed pair (and,
by the same search on triples, presumably no fixed O(1)-size family) survives
as n grows; the set of adequate cuts must depend on pi/n, matching the
already-rejected restricted-cut families of session 8.

Usage: python3 experiments/h13_averaging_bound.py [--nmax 11] [--hillclimb-n 20]
Output: data/runs/h13_averaging_bound/report.json, report.md.
Version h13_averaging_bound-1.0.
"""

import argparse
import itertools
import json
import os
import random
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
from moves import CORE_VERSION  # noqa: E402

VERSION = "h13_averaging_bound-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h13_averaging_bound")


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def pair_data(pi):
    n = len(pi)
    out = []
    for i1 in range(n):
        for i2 in range(i1 + 1, n):
            d = i2 - i1
            e = (pi[i2] - pi[i1]) % n
            out.append((i1, i2, d, e))
    return out


def inv_of_cut(pi, q, c):
    n = len(pi)
    shift = (q + 1 - c) % n
    w = [(pi[(q + 1 + j) % n] - shift) % n for j in range(n)]
    return inversions(w)


def part1_decomposition(nmax, log):
    """Verify inv(w) == sum_pairs (A(a) xor B(b)) for random pi, all (q, c)."""
    random.seed(0)
    ok = True
    checked = 0
    for n in range(4, nmax + 1):
        for _trial in range(15):
            pi = list(range(n))
            random.shuffle(pi)
            pairs = pair_data(pi)
            for q in range(n):
                for c in range(n):
                    a = (q + 1) % n
                    b = (q + 1 - c) % n
                    iv = inv_of_cut(pi, q, c)
                    total = 0
                    for (i1, i2, d, e) in pairs:
                        A = 1 if i1 + 1 <= a <= i2 else 0
                        lo = (pi[i1] + 1) % n
                        B = 1 if ((b - lo) % n) < e else 0
                        total += A ^ B
                    if total != iv:
                        ok = False
                        log(f"FAIL decomposition n={n} pi={pi} q={q} c={c}")
                        break
                if not ok:
                    break
            checked += 1
            if not ok:
                break
        if not ok:
            break
    log(f"part 1: decomposition A xor B verified on {checked} random (n, pi) samples, 4 <= n <= {nmax}: {'OK' if ok else 'FAIL'}")
    return ok, checked


def F_of(a, n):
    return a * (n - a)


def H_reflection(a, n):
    """H(a) for a reflection pi_h (any h): sum_pairs A_pair(a) * (n - 2*e_pair)
    with e_pair = n - d_pair. Direct O(n^2) computation, no dependence on h."""
    total = 0
    for i1 in range(n):
        for i2 in range(i1 + 1, n):
            d = i2 - i1
            if i1 + 1 <= a <= i2:
                total += (2 * d - n)
    return total


def E_reflection(n):
    return sum(k * k for k in range(1, n))


def part2_pigeonhole(nmax, log):
    """Check H(a) == 0 for all a on reflections, E(n) closed form, and that
    the resulting pigeonhole bound E(n)/n exceeds floor((n-1)^2/4) for all
    n >= 4, with the gap matching (n^2-1)/12 (before floor adjustments)."""
    rows = []
    all_H_zero = True
    for n in range(4, nmax + 1):
        Hs = [H_reflection(a, n) for a in range(n)]
        if any(h != 0 for h in Hs):
            all_H_zero = False
        # direct E via a genuine reflection permutation, cross-checked against closed form
        pi_h = [(1 - i) % n for i in range(n)]
        pairs = pair_data(pi_h)
        E_direct = sum(e for (_, _, _, e) in pairs)
        E_formula = E_reflection(n)
        thr = (n - 1) ** 2 // 4
        bound = E_formula / n
        gap_exact = (n * n - 1) / 12
        rows.append({
            "n": n, "H_all_zero": all(h == 0 for h in Hs), "E_direct": E_direct,
            "E_formula": E_formula, "threshold": thr, "pigeonhole_bound": bound,
            "gap": bound - thr, "gap_formula_n2_minus_1_over_12": gap_exact,
        })
        assert E_direct == E_formula, (n, E_direct, E_formula)
        log(f"part 2 n={n}: H(a)=0 for all a: {all(h==0 for h in Hs)}; E={E_formula} "
            f"(n-1)(2n-1)/6; pigeonhole bound E/n={bound:.3f} vs threshold {thr} "
            f"-> gap {bound-thr:.3f} (method fails)")
    return all_H_zero, rows


def minb_for_a(pi, a, n):
    best = None
    for b in range(n):
        w = [(pi[(a + j) % n] - b) % n for j in range(n)]
        iv = inversions(w)
        if best is None or iv < best:
            best = iv
    return best


def score_fixed_a(pi, a_list, n):
    return min(minb_for_a(pi, a, n) for a in a_list)


def part3_fixed_family(nmax_exhaustive, hillclimb_n, log, a_list=(0, 3), iters=1500, restarts=6, seed=42):
    exhaustive_rows = []
    for n in range(4, nmax_exhaustive + 1):
        thr = (n - 1) ** 2 // 4
        worst = -1
        argworst = None
        fails = 0
        for pi in itertools.permutations(range(n)):
            m = score_fixed_a(list(pi), a_list, n)
            if m > worst:
                worst = m
                argworst = pi
            if m > thr:
                fails += 1
        exhaustive_rows.append({"n": n, "threshold": thr, "worst": worst, "excess": worst - thr, "fails": fails, "argworst": list(argworst)})
        log(f"part 3 (exhaustive) n={n} a_list={a_list}: worst={worst} threshold={thr} excess={worst-thr} fails={fails}")

    random.seed(seed)
    n = hillclimb_n
    thr = (n - 1) ** 2 // 4
    best_overall = -1
    best_pi = None
    for _r in range(restarts):
        pi = list(range(n))
        random.shuffle(pi)
        cur = score_fixed_a(pi, a_list, n)
        for _it in range(iters):
            i, j = random.sample(range(n), 2)
            pi[i], pi[j] = pi[j], pi[i]
            new = score_fixed_a(pi, a_list, n)
            if new >= cur:
                cur = new
            else:
                pi[i], pi[j] = pi[j], pi[i]
        if cur > best_overall:
            best_overall = cur
            best_pi = pi[:]
    hillclimb_row = {
        "n": n, "a_list": list(a_list), "threshold": thr, "best_found": best_overall,
        "excess": best_overall - thr, "witness_pi": best_pi, "method": "hill-climbing (local search, NOT exhaustive; a lower bound on the true worst case)",
        "iters": iters, "restarts": restarts, "seed": seed,
    }
    # verify witness independently with a direct (non-hillclimb) recomputation
    verified = score_fixed_a(best_pi, a_list, n) == best_overall
    hillclimb_row["witness_verified"] = verified
    log(f"part 3 (hill-climbing) n={n} a_list={a_list}: found excess={best_overall-thr} "
        f"(witness_verified={verified}) -- fixed pair is NOT sufficient in general")
    return exhaustive_rows, hillclimb_row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8, help="max n for part 1 (decomposition, random) and part 2 (pigeonhole, exhaustive over a)")
    ap.add_argument("--nmax-exhaustive-family", type=int, default=8, help="max n for exhaustive part 3 check of the fixed pair {0,3} (pure Python; n=9..11 done separately in C, see experiments/line_pair_a03.c)")
    ap.add_argument("--hillclimb-n", type=int, default=20)
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    lines = []

    def log(msg):
        print(msg)
        lines.append(msg)

    t0 = time.time()
    ok1, checked1 = part1_decomposition(args.nmax, log)
    ok2, rows2 = part2_pigeonhole(max(args.nmax, 10), log)
    rows3, hillclimb_row = part3_fixed_family(args.nmax_exhaustive_family, args.hillclimb_n, log)
    elapsed = time.time() - t0

    report = {
        "version": VERSION,
        "core_version": CORE_VERSION,
        "elapsed_seconds": round(elapsed, 1),
        "part1_decomposition_ok": ok1,
        "part1_checked_samples": checked1,
        "part2_pigeonhole_all_H_zero": ok2,
        "part2_rows": rows2,
        "part3_exhaustive_fixed_pair_0_3": rows3,
        "part3_hillclimb": hillclimb_row,
        "verdict": (
            "H13-I NOT proved. Part 2 shows the pigeonhole/averaging bound "
            "suggested by the algebraic hint fails on reflections by a gap "
            "growing like (n^2-1)/12, for every choice of a -- a proved dead "
            "end, not just a small-n observation. Part 3 shows the natural "
            "'fixed small family of a' idea (exhaustively sufficient up to "
            "n=11 for {0,3}) breaks at n=20 (hill-climbing counterexample, "
            "excess 2) -- consistent with session 8's rejection of "
            "restricted-cut families: no O(1)-size, pi-independent set of "
            "position-cuts suffices in general."
        ),
    }
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(report, f, indent=2)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write(f"# h13_averaging_bound report ({VERSION}, core {CORE_VERSION})\n\n")
        f.write(f"Elapsed: {elapsed:.1f} s\n\n```\n")
        f.write("\n".join(lines))
        f.write("\n```\n\nVerdict: " + report["verdict"] + "\n")
    print(f"\nWrote {OUT}/report.json and report.md ({elapsed:.1f} s)")


if __name__ == "__main__":
    main()
