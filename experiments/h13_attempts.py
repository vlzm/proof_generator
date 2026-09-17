"""H13-I, session 9: two attempted routes to a proof, both checked and rejected.

Route 1 (single-axis rotation). Restrict the double cut (q, c) to one of the
two n-sized "diagonals" q fixed / c fixed (only rotate values, or only rotate
positions). Exhaustively, for 4 <= n <= AMAX, compute
    max_pi min_c inv(pi read in original position order, values rotated by c)
and the symmetric max_pi min_q inv(pi rotated in position, values as-is), and
compare to floor((n-1)^2/4). Both routes fail starting at n = 7 (the same
permutation (0,5,4,3,2,1,6) that already refutes Lemma C in
docs/notes/h13_line_model.md Sec. 4 -- unsurprising, since fixing one cut is
a strictly smaller search than the n^2 double cut).

Route 2 (second-moment / Cantelli bound on the C37 distribution). For a
permutation pi, let X be inv(w_{q,c}) as (q, c) ranges uniformly over the n^2
cuts (mean mu = S(pi)/n^2 from C37, max M, variance Var). The elementary
inequality Var(X) <= (M - mu)(mu - m) for any bounded random variable
(m = min X) gives m <= mu - Var(X) / (M - mu) whenever M > mu. This route
computes that guaranteed upper bound ("Cantelli estimate") on rev_n for
4 <= n <= RMAX and compares it to floor((n-1)^2/4) and to the true I(rev_n).
The estimate exceeds the target bound already at n = 4 and the gap grows with
n, so this route also fails (it is a valid inequality, just not tight enough
here -- it does not by itself say the claim is false, since the true minimum
is checked separately and matches floor((n-1)^2/4) throughout).

Neither route proves or refutes H13-I; see docs/notes/h13_line_model.md Sec. 7.
Usage: python3 experiments/h13_attempts.py [--amax 8] [--rmax 12]
Output: data/runs/h13_attempts/report.json, report.md.
"""

import argparse
import itertools
import json
import os
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
from moves import CORE_VERSION  # noqa: E402

VERSION = "h13_attempts-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h13_attempts")


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def value_rotation_only_min(pi):
    n = len(pi)
    return min(inversions([(pi[i] - c) % n for i in range(n)]) for c in range(n))


def position_rotation_only_min(pi):
    n = len(pi)
    return min(inversions([pi[(q + 1 + j) % n] for j in range(n)]) for q in range(n))


def all_inv_list(pi):
    n = len(pi)
    vals = []
    for q in range(n):
        for c in range(n):
            shift = (q + 1 - c) % n
            w = [(pi[(q + 1 + j) % n] - shift) % n for j in range(n)]
            vals.append(inversions(w))
    return vals


def route1(amax, log):
    rows = []
    ok_route = True  # "ok" here means the restricted route WOULD have sufficed (expected False)
    for n in range(4, amax + 1):
        bound = (n - 1) ** 2 // 4
        worst_v = 0
        worst_v_pi = None
        worst_p = 0
        worst_p_pi = None
        for pi in itertools.permutations(range(n)):
            v = value_rotation_only_min(pi)
            if v > worst_v:
                worst_v, worst_v_pi = v, pi
            p = position_rotation_only_min(pi)
            if p > worst_p:
                worst_p, worst_p_pi = p, pi
        fails = worst_v > bound or worst_p > bound
        ok_route = ok_route and not fails
        rows.append({"n": n, "bound": bound, "max_value_rot_only": worst_v, "argmax_value": worst_v_pi,
                     "max_position_rot_only": worst_p, "argmax_position": worst_p_pi})
        log(f"route1 n={n}: bound={bound}, max over pi of value-rotation-only min = {worst_v} "
            f"(pi={worst_v_pi}), max over pi of position-rotation-only min = {worst_p} (pi={worst_p_pi}) "
            f"-> {'both within bound' if not fails else 'EXCEEDS bound (route fails)'}")
    return ok_route, rows


def route2(rmax, log):
    rows = []
    for n in range(4, rmax + 1):
        pi = tuple(range(n - 1, -1, -1))  # rev_n, a reflection: known I(rev_n) = floor((n-1)^2/4)
        vals = all_inv_list(pi)
        m, M, N = min(vals), max(vals), len(vals)
        mu = sum(vals) / N
        var = sum((v - mu) ** 2 for v in vals) / N
        est = mu - var / (M - mu) if M > mu else mu
        bound = (n - 1) ** 2 // 4
        rows.append({"n": n, "bound": bound, "true_I_rev_n": m, "mean": round(mu, 3),
                     "max": M, "cantelli_estimate": round(est, 3)})
        log(f"route2 n={n}: bound={bound}, true I(rev_n)={m}, mean={mu:.2f}, max={M}, "
            f"Cantelli guaranteed upper bound={est:.2f} -> "
            f"{'would prove it' if est <= bound + 1e-9 else 'NOT tight enough'}")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amax", type=int, default=8)
    ap.add_argument("--rmax", type=int, default=12)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    lines = []

    def log(msg):
        print(msg, flush=True)
        lines.append(msg)

    log(f"== {VERSION} {CORE_VERSION} args={vars(args)}")
    t0 = time.time()
    ok1, r1 = route1(args.amax, log)
    r2 = route2(args.rmax, log)
    log(f"route1 sufficed for all n checked: {ok1} (expected False -- both single-axis "
        f"restrictions fail already at n=7)")
    log(f"total time: {time.time() - t0:.1f} s")

    report = {"version": VERSION, "core": CORE_VERSION, "args": vars(args),
              "route1_rotation_only": r1, "route1_would_suffice": ok1,
              "route2_cantelli": r2}
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(report, f, indent=2)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# h13_attempts report (session 9)\n\n```\n" + "\n".join(lines) + "\n```\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
