"""H13-I, session 9: torus/quadrant reformulation and three rejected proof attempts.

H13-I claims I(pi) = min_{a,b in Z_n} inv(w(a,b)) <= floor((n-1)^2/4) for every
permutation pi of Z_n, where w(a,b)_j = (pi((a+j) mod n) - b) mod n, j = 0..n-1
(a = q+1 is the position-cut rotation, b = q+1-c is the value-cut rotation;
this is the same I(pi) as C33/line_profile, re-derived here with independent,
elementary code -- no dependency on oracle/moves.py or exact/bfs.py).

Part 1 (reformulation, PROVED): concordant(a,b) := C(n,2) - inv(w(a,b)) counts
pairs {i, i'} of positions whose cyclic order read from a agrees with the
cyclic order of their values read from b ("same-quadrant pairs" of the task's
torus framing).  H13-I <=> max_{a,b} concordant(a,b) >= floor(n^2/4) (via the
identity floor(n^2/4) + floor((n-1)^2/4) = n(n-1)/2 = C(n,2), C12).  For an
ordered pair (i, i') with position gap g = (i'-i) mod n in [1,n-1] and value
gap h = (pi(i')-pi(i)) mod n in [1,n-1], exactly (n-g) values of a put i before
i' (and g values put i' before i); symmetrically for b and h.  Hence the
number of (a,b) in Z_n x Z_n making the pair {i,i'} concordant is exactly
  K(g,h) = (n-g)(n-h) + g*h.
Summing over the n(n-1)/2 unordered pairs and dividing by n^2 recovers the
"averaging over all n^2 cuts" bound already rejected in session 8 (insufficient
in both directions: h13_line_model.md sec 6).  check_reformulation() verifies
sum_{a,b} concordant(a,b) == sum_{pairs} K(g,h) exhaustively.

Part 2 (second-moment / Cauchy-Schwarz bound, REFUTED as a proof technique):
max_{a,b} concordant(a,b) >= S2/S where S = sum concordant(a,b), S2 = sum
concordant(a,b)^2 over the n^2 grid (since max*S >= sum x_i^2).  This uses
more information than the plain mean and could in principle beat it, but it
already fails at the exact equality case: n=4, pi=(0,3,2,1) has max = 4 =
floor(n^2/4) (tight) while S2/S = 17/5 = 3.4 < 4.  check_second_moment() finds
the worst case for each n.

Part 3 (coupled/joint local optimality, REFUTED as sufficient): call (a,b) a
"joint fixed point" if b globally minimises inv(w(a,.)) over its whole row AND
a globally minimises inv(w(.,b)) over its whole column (this is strictly
stronger than Lemma C in h13_line_model.md, which only checks unit and
k-step shifts along the two elementary directions (a+k,b+k) and (a,b-k)
separately, never a genuinely joint move).  The global minimiser is always
such a fixed point, so existence is trivial; the question is whether *every*
fixed point is good.  It is not: n=7, pi=(0,5,4,3,2,1,6) (the Lemma C
counterexample of h13_line_model.md sec 4) has fixed point (a,b)=(0,0) with
inv=10 > floor(6^2/4)=9, coexisting with good fixed points such as (1,2) with
inv=8.  check_joint_fixed_points() finds all such bad-fixed-point permutations
exhaustively for small n.

Part 4 (321-avoidance / Mantel bound, REFUTED): floor((n-1)^2/4) is the Turan
(Mantel) bound for triangle-free graphs on n-1 vertices, suggesting the
inversion graph of an optimal cut might always be triangle-free (i.e. w(a*,b*)
avoids the pattern 321, splitting into <=2 increasing runs).  False already at
n=5: pi=(0,1,4,3,2) has best inv=3, and every optimal cut's line contains a
321 pattern.  check_321_avoidance() counts how often this happens.

None of parts 2-4 give a proof; see docs/notes/h13_i_session9.md for the
verdict.  Usage: python3 experiments/h13_i_quadrant.py --nmax 8
Output: data/runs/h13_i_session9/report.json report.md.  Version h13_i_quadrant-1.0.
"""

import argparse
import itertools
import json
import os
import time

VERSION = "h13_i_quadrant-1.0"
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "data", "runs", "h13_i_session9")


def inv_count(w):
    n = len(w)
    c = 0
    for i in range(n):
        for j in range(i + 1, n):
            if w[i] > w[j]:
                c += 1
    return c


def line(pi, n, a, b):
    return tuple((pi[(a + j) % n] - b) % n for j in range(n))


def inv_ab(pi, n, a, b):
    return inv_count(line(pi, n, a, b))


def full_grid(pi, n):
    """inv(a, b) for all a, b in Z_n; independent elementary O(n^4) code."""
    return [[inv_ab(pi, n, a, b) for b in range(n)] for a in range(n)]


def check_reformulation(nmax):
    """Verify sum_{a,b} concordant(a,b) == sum_{pairs} K(g,h) exhaustively."""
    results = []
    for n in range(4, nmax + 1):
        C = n * (n - 1) // 2
        ok = True
        checked = 0
        for pi in itertools.permutations(range(n)):
            grid = full_grid(pi, n)
            total_direct = sum(C - grid[a][b] for a in range(n) for b in range(n))
            total_formula = 0
            for i in range(n):
                for ip in range(i + 1, n):
                    g = (ip - i) % n
                    h = (pi[ip] - pi[i]) % n
                    total_formula += (n - g) * (n - h) + g * h
            if total_direct != total_formula:
                ok = False
            checked += 1
        results.append({"n": n, "perms_checked": checked, "identity_holds": ok})
        print(f"  n={n}: identity K(g,h) verified on all {checked} perms: {ok}")
    return results


def check_second_moment(nmax):
    """Find, for each n, the permutation minimising S2/S - floor(n^2/4)."""
    results = []
    for n in range(4, nmax + 1):
        target = n * n // 4
        C = n * (n - 1) // 2
        worst_gap = None
        worst_pi = None
        worst_bound = None
        worst_max = None
        for pi in itertools.permutations(range(n)):
            grid = full_grid(pi, n)
            vals = [C - grid[a][b] for a in range(n) for b in range(n)]
            S = sum(vals)
            S2 = sum(v * v for v in vals)
            bound = S2 / S if S else 0.0
            gap = bound - target
            if worst_gap is None or gap < worst_gap:
                worst_gap = gap
                worst_pi = pi
                worst_bound = bound
                worst_max = max(vals)
        fails = worst_gap < -1e-9
        results.append({
            "n": n, "target": target, "worst_pi": list(worst_pi),
            "S2_over_S": worst_bound, "actual_max": worst_max,
            "second_moment_bound_holds": not fails,
        })
        print(f"  n={n}: worst pi={worst_pi} S2/S={worst_bound:.3f} target={target} "
              f"actual_max={worst_max} bound_holds={not fails}")
    return results


def best_b_for_a(pi, n, a):
    vals = [inv_ab(pi, n, a, b) for b in range(n)]
    m = min(vals)
    return m


def check_joint_fixed_points(nmax):
    """Count, for each n, permutations whose worst joint fixed point exceeds
    floor((n-1)^2/4); a fixed point is globally optimal on its whole row AND
    its whole column (Part 3 above)."""
    results = []
    for n in range(4, nmax + 1):
        target = (n - 1) * (n - 1) // 4
        bad = []
        total = 0
        for pi in itertools.permutations(range(n)):
            total += 1
            grid = full_grid(pi, n)
            row_min = [min(grid[a]) for a in range(n)]
            col_min = [min(grid[a][b] for a in range(n)) for b in range(n)]
            fps = [(a, b, grid[a][b]) for a in range(n) for b in range(n)
                   if grid[a][b] == row_min[a] and grid[a][b] == col_min[b]]
            worst_fp = max(v for _, _, v in fps)
            if worst_fp > target:
                bad.append((pi, worst_fp, fps))
        results.append({
            "n": n, "target": target, "perms_total": total,
            "perms_with_bad_fixed_point": len(bad),
            "example": ({
                "pi": list(bad[0][0]), "worst_fixed_point_inv": bad[0][1],
                "all_fixed_points": bad[0][2],
            } if bad else None),
        })
        print(f"  n={n}: {len(bad)}/{total} perms have a bad joint fixed point (target {target})")
    return results


def has_321(w):
    """True iff w contains a decreasing subsequence of length 3 (patience sort)."""
    import bisect
    tails = []
    for x in w:
        y = -x
        idx = bisect.bisect_left(tails, y)
        if idx == len(tails):
            tails.append(y)
        else:
            tails[idx] = y
        if len(tails) >= 3:
            return True
    return False


def check_321_avoidance(nmax):
    """Count permutations for which NO optimal cut's line is 321-avoiding."""
    results = []
    for n in range(4, nmax + 1):
        total = 0
        no_avoiding_opt = 0
        example = None
        for pi in itertools.permutations(range(n)):
            total += 1
            best = None
            best_lines = []
            for a in range(n):
                for b in range(n):
                    w = line(pi, n, a, b)
                    iv = inv_count(w)
                    if best is None or iv < best:
                        best = iv
                        best_lines = [w]
                    elif iv == best:
                        best_lines.append(w)
            if all(has_321(w) for w in best_lines):
                no_avoiding_opt += 1
                if example is None:
                    example = {"pi": list(pi), "best_inv": best, "an_optimal_line": list(best_lines[0])}
        results.append({
            "n": n, "perms_total": total,
            "perms_where_no_optimal_cut_avoids_321": no_avoiding_opt,
            "example": example,
        })
        print(f"  n={n}: {no_avoiding_opt}/{total} perms have no 321-avoiding optimal cut")
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=7, help="exhaustive up to this n (cost is n! * n^4)")
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()

    print("Part 1: reformulation identity sum concordant(a,b) = sum K(g,h)")
    part1 = check_reformulation(args.nmax)
    print("Part 2: second-moment (Cauchy-Schwarz) bound on max concordant(a,b)")
    part2 = check_second_moment(args.nmax)
    print("Part 3: joint (row+column) fixed points of inv(a,b)")
    part3 = check_joint_fixed_points(args.nmax)
    print("Part 4: 321-avoidance of optimal cuts")
    part4 = check_321_avoidance(args.nmax)

    elapsed = time.time() - t0
    report = {
        "version": VERSION, "nmax": args.nmax, "elapsed_s": elapsed,
        "part1_reformulation_identity": part1,
        "part2_second_moment_bound": part2,
        "part3_joint_fixed_points": part3,
        "part4_321_avoidance": part4,
    }
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(report, f, indent=2)

    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write(f"# H13-I session 9: quadrant reformulation and rejected attempts ({VERSION})\n\n")
        f.write(f"nmax = {args.nmax}, elapsed = {elapsed:.1f} s\n\n")
        f.write("## Part 1: sum_{a,b} concordant(a,b) = sum_{pairs} K(g,h) -- identity, PROVED\n\n")
        for r in part1:
            f.write(f"- n={r['n']}: {r['perms_checked']} perms, identity holds: {r['identity_holds']}\n")
        f.write("\n## Part 2: second-moment bound S2/S -- REFUTED (fails already at the equality case)\n\n")
        for r in part2:
            f.write(f"- n={r['n']}: target={r['target']}, worst pi={r['worst_pi']}, "
                    f"S2/S={r['S2_over_S']:.3f}, actual max={r['actual_max']}, "
                    f"bound holds: {r['second_moment_bound_holds']}\n")
        f.write("\n## Part 3: joint (row+column) fixed points -- REFUTED as sufficient\n\n")
        for r in part3:
            f.write(f"- n={r['n']}: {r['perms_with_bad_fixed_point']}/{r['perms_total']} perms have a "
                    f"joint fixed point exceeding floor((n-1)^2/4)={r['target']}\n")
            if r["example"]:
                f.write(f"  example: pi={r['example']['pi']}, "
                        f"worst fixed point inv={r['example']['worst_fixed_point_inv']}, "
                        f"all fixed points={r['example']['all_fixed_points']}\n")
        f.write("\n## Part 4: 321-avoidance of optimal cuts -- REFUTED\n\n")
        for r in part4:
            f.write(f"- n={r['n']}: {r['perms_where_no_optimal_cut_avoids_321']}/{r['perms_total']} perms "
                    f"have no 321-avoiding optimal cut\n")
            if r["example"]:
                f.write(f"  example: pi={r['example']['pi']}, best inv={r['example']['best_inv']}, "
                        f"an optimal line={r['example']['an_optimal_line']}\n")

    print(f"Done in {elapsed:.1f}s. Report: {OUT}/report.md")


if __name__ == "__main__":
    main()
