"""H13-I candidate strategies (session 9): does a cheaper rule than the full
n^2-cut minimum already witness I(pi) <= floor((n-1)^2/4)?

Background: docs/notes/h13_line_model.md #6. I(pi) = min over n^2 cuts (q,s)
of the number of inversions of the relabelled line (see that note's #0 for
the exact definition of a cut and the line w). This script tests three
candidate ways to find a good cut cheaply and the Mantel/Turan reformulation
of the target bound; see docs/notes/h13i_session9.md for the write-up.

Definitions (all exhaustive over all n^2 cuts unless noted):
  I_full(pi)      -- true I(pi), min over all (q, s).
  minq_avgs(pi)    -- min over q of (average over s of inv(q, s)), i.e. the
                      best q under the *expected* inversions over s; a lower
                      bound witness only if it is itself <= floor((n-1)^2/4)
                      (pigeonhole: some s attains at most the average).
  anchor_min(pi)   -- min over the n cuts (q, s) = (i0, pi(i0)) that anchor
                      the origin at a data point of pi itself (induction-style
                      reduction to n-1 elements, PLAN #8 / notes #1.7).
  best_align(pi)   -- max over the n position-arcs A of size ceil(n/2) and
                      the n value-arcs B of size ceil(n/2) (n^2 pairs) of
                      x = |A cap pi^{-1}(B)|; the guaranteed-concordant count
                      from a perfect (A, B) alignment is x * (n - 2*ceil(n/2) + x)
                      (see docs/notes/h13i_session9.md); the pair is "perfect"
                      when x = ceil(n/2).

Usage: python3 experiments/h13i_candidates.py --nmax 8
Output: data/runs/h13i_candidates/report.md, report.json. Version 1.0.
"""

import argparse
import itertools
import json
import os
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def inv_count(seq):
    n = len(seq)
    c = 0
    for i in range(n):
        si = seq[i]
        for j in range(i + 1, n):
            if si > seq[j]:
                c += 1
    return c


def inv_at_cut(pi, q, s, n):
    order = [pi[(q + 1 + j) % n] for j in range(n)]
    w = [(v - s) % n for v in order]
    return inv_count(w)


def I_full(pi, n):
    best = None
    for q in range(n):
        for s in range(n):
            ic = inv_at_cut(pi, q, s, n)
            if best is None or ic < best:
                best = ic
    return best


def F_of_q(pi, q, n):
    order = [pi[(q + 1 + j) % n] for j in range(n)]
    total = 0
    for j in range(n):
        for k in range(j + 1, n):
            total += (order[k] - order[j]) % n
    return total


def minq_avgs(pi, n):
    # min over q of sum_s inv(q, s) / n; report the (unrounded) sum, compare
    # against n * bound.
    return min(F_of_q(pi, q, n) for q in range(n))


def anchor_min(pi, n):
    best = None
    for q in range(n):
        s = pi[q]
        ic = inv_at_cut(pi, q, s, n)
        if best is None or ic < best:
            best = ic
    return best


def best_align(pi, n):
    a = (n + 1) // 2
    inv_pi = [0] * n
    for i in range(n):
        inv_pi[pi[i]] = i
    best_x = -1
    for qpos in range(n):
        pos_arc = set((qpos + j) % n for j in range(a))
        for qval in range(n):
            val_arc = set((qval + j) % n for j in range(a))
            x = sum(1 for i in pos_arc if pi[i] in val_arc)
            if x > best_x:
                best_x = x
    return best_x, a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=4)
    ap.add_argument("--nmax", type=int, default=8)
    args = ap.parse_args()

    out_dir = os.path.join(ROOT, "data", "runs", "h13i_candidates")
    os.makedirs(out_dir, exist_ok=True)

    results = {}
    for n in range(args.nmin, args.nmax + 1):
        t0 = time.time()
        bound = (n - 1) ** 2 // 4
        a = (n + 1) // 2
        max_I = 0
        max_minq_avgs = -1
        worst_minq_avgs_pi = None
        max_anchor = 0
        worst_anchor_pi = None
        min_best_align_x = n + 1
        worst_align_pi = None
        n_perfect_align = 0
        checked = 0
        for pi in itertools.permutations(range(n)):
            checked += 1
            I = I_full(pi, n)
            if I > max_I:
                max_I = I
            mqa = minq_avgs(pi, n)
            if mqa > max_minq_avgs:
                max_minq_avgs = mqa
                worst_minq_avgs_pi = pi
            am = anchor_min(pi, n)
            if am > max_anchor:
                max_anchor = am
                worst_anchor_pi = pi
            x, _ = best_align(pi, n)
            if x < min_best_align_x:
                min_best_align_x = x
                worst_align_pi = pi
            if x == a:
                n_perfect_align += 1
        dt = time.time() - t0
        results[n] = dict(
            bound=bound,
            n_bound=n * bound,
            max_I=max_I,
            max_minq_avgs=max_minq_avgs,
            worst_minq_avgs_pi=list(worst_minq_avgs_pi),
            max_anchor=max_anchor,
            worst_anchor_pi=list(worst_anchor_pi),
            a=a,
            min_best_align_x=min_best_align_x,
            worst_align_pi=list(worst_align_pi),
            n_perfect_align=n_perfect_align,
            n_perms=checked,
            seconds=dt,
        )
        print(
            f"n={n} bound={bound} max_I={max_I} "
            f"max_minq_avgs={max_minq_avgs} (n*bound={n*bound}) "
            f"max_anchor={max_anchor} (bound={bound}) "
            f"min_best_align_x={min_best_align_x} (need {a}) "
            f"perfect_align={n_perfect_align}/{checked} "
            f"[{dt:.1f}s]"
        )

    with open(os.path.join(out_dir, "report.json"), "w") as f:
        json.dump(results, f, indent=2)

    with open(os.path.join(out_dir, "report.md"), "w") as f:
        f.write("# H13-I candidate strategies -- report\n\n")
        f.write(
            "Exhaustive over all pi, n = %d..%d. Version h13i_candidates-1.0.\n\n"
            % (args.nmin, args.nmax)
        )
        f.write(
            "| n | bound floor((n-1)^2/4) | max I(pi) | max minq_avgs "
            "(vs n*bound) | max anchor_min (vs bound) | min best_align x "
            "(need ceil(n/2)) | perfect alignments / perms | time |\n"
        )
        f.write("|---|---|---|---|---|---|---|---|\n")
        for n, r in results.items():
            f.write(
                "| %d | %d | %d | %d (%d) | %d (%d) | %d (%d) | %d/%d | %.1fs |\n"
                % (
                    n,
                    r["bound"],
                    r["max_I"],
                    r["max_minq_avgs"],
                    r["n_bound"],
                    r["max_anchor"],
                    r["bound"],
                    r["min_best_align_x"],
                    r["a"],
                    r["n_perfect_align"],
                    r["n_perms"],
                    r["seconds"],
                )
            )
        f.write("\nWorst cases (per n):\n\n")
        for n, r in results.items():
            f.write(
                "- n=%d: minq_avgs worst pi = %s; anchor worst pi = %s; "
                "align worst pi = %s\n"
                % (
                    n,
                    r["worst_minq_avgs_pi"],
                    r["worst_anchor_pi"],
                    r["worst_align_pi"],
                )
            )


if __name__ == "__main__":
    main()
