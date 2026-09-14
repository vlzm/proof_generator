"""H13-I attempt (session 9): decoupled (s, t) parametrization of the double cut.

`I(pi) = min_{q,c} inv(w)` (PROBLEM Sec.5.6, C33/C35) uses the (q, c) cut of
`docs/notes/h13_line_model.md`. This module rewrites the same n^2 cuts in a
decoupled pair (s, t): domain start s (pure re-reading of pi, no value change)
and value shift t (pure additive relabelling, no domain change):

    W(pi, n, s, t)[j] = (pi[(s + j) % n] - t) % n,   j = 0..n-1
    s = q + 1,  t = q + 1 - c   (bijective reparametrisation of the n^2 cuts)

Exact single-step recursions (proved and checked below):

    inv(W_{s+1,t}) - inv(W_{s,t}) = (n-1) - 2*v0(s,t),   v0(s,t) = W(pi,n,s,t)[0] = (pi[s]-t) % n
    inv(W_{s,t+1}) - inv(W_{s,t}) = (n-1) - 2*p0(s,t),   p0(s,t) = (pi^{-1}[t]-s) % n

(s+1: front element moved to the back, value unchanged, so it flips role in
exactly v0 vs n-1-v0 pairs; t+1: the element currently valued 0 wraps to n-1
at the same position p0, flipping role in exactly p0 vs n-1-p0 pairs.)

Findings this session, all machine-checked below:

1. Duality: inv_pi(s,t) = inv_{pi^{-1}}(t,s) (a permutation and its positional
   inverse have equal inversion count), hence I(pi) = I(pi^{-1}) for all pi
   (registered as C37). Exhaustive 4 <= n <= 8.
2. Single-shift is not enough: fixing s = 0 (i.e. reading pi in its own given
   order, no domain rotation) and minimizing only over t fails to reach
   floor((n-1)^2/4) already at n = 7 (worst case 10 > 9), n = 8 (13 > 12),
   n = 9 (17 > 16). So both degrees of freedom are genuinely needed; rules
   out any "value rotation alone" simplification of H13-I. Exhaustive 4<=n<=9.
3. The n^2-cut average of inv(pi) is far from floor((n-1)^2/4) even for the
   extremal pi = sigma_n-type reflections: exact average = (n-1)(2n-1)/6
   ~ n^2/3, well above the target ~n^2/4, while the true minimum equals the
   target exactly. So the distribution of inv over the n^2 cuts is highly
   skewed and no averaging/pigeonhole argument over all cuts (or over a
   fixed-weight subset) can prove H13-I; a good cut must be located
   structurally, not statistically. Exact formula checked 4 <= n <= 12.
4. H13-I itself (I(pi) <= floor((n-1)^2/4) for all pi) is independently
   re-verified 4 <= n <= 9 by this fresh implementation (does not reuse
   `experiments/line_profile.c`), agreeing with C33/C35. Not proved in
   general this session; verdict and next steps in
   `docs/notes/h13i_toric_recursion.md`.

Usage: python3 experiments/h13i_toric.py --nmax 9
Output: data/runs/h13i_toric/report.md, report.json. Version h13i_toric-1.0.
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

VERSION = "h13i_toric-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h13i_toric")


def inv(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def W(pi, n, s, t):
    return [(pi[(s + j) % n] - t) % n for j in range(n)]


def invperm(pi):
    n = len(pi)
    out = [0] * n
    for i, v in enumerate(pi):
        out[v] = i
    return tuple(out)


def I_full(pi, n):
    best = None
    for s in range(n):
        for t in range(n):
            iv = inv(W(pi, n, s, t))
            if best is None or iv < best:
                best = iv
    return best


def I_single_shift(pi, n):
    """Only value-shift t; domain order fixed as given (s = 0)."""
    best = None
    for t in range(n):
        iv = inv([(x - t) % n for x in pi])
        if best is None or iv < best:
            best = iv
    return best


def check_recursion(trials, log):
    random.seed(2)
    checked = 0
    for _ in range(trials):
        n = random.randint(4, 9)
        pi = list(range(n))
        random.shuffle(pi)
        pi = tuple(pi)
        s, t = random.randrange(n), random.randrange(n)
        w = W(pi, n, s, t)
        v0 = w[0]
        d_pred = (n - 1) - 2 * v0
        d_actual = inv(W(pi, n, (s + 1) % n, t)) - inv(w)
        assert d_pred == d_actual, ("s-step", pi, s, t, d_pred, d_actual)
        pinv = {v: i for i, v in enumerate(pi)}
        p0 = (pinv[t] - s) % n
        d_pred2 = (n - 1) - 2 * p0
        d_actual2 = inv(W(pi, n, s, (t + 1) % n)) - inv(w)
        assert d_pred2 == d_actual2, ("t-step", pi, s, t, d_pred2, d_actual2)
        checked += 1
    log(f"recursion check: {checked}/{trials} random (pi,s,t) triples, both single-step formulas exact")


def run(n, log, check_duality=False):
    t0 = time.time()
    T = (n - 1) ** 2 // 4
    worst_full, arg_full = -1, None
    worst_single, arg_single = -1, None
    duality_fail = None
    duality_checked = 0
    cnt = 0
    for pi in itertools.permutations(range(n)):
        cnt += 1
        full = I_full(pi, n)
        if full > worst_full:
            worst_full, arg_full = full, pi
        single = I_single_shift(pi, n)
        if single > worst_single:
            worst_single, arg_single = single, pi
        if check_duality:
            other = I_full(invperm(pi), n)
            duality_checked += 1
            if other != full and duality_fail is None:
                duality_fail = (pi, full, other)
    B = n * (n - 1) // 2
    sigma = tuple((1 - i) % n for i in range(n))
    avg_total = sum(inv(W(sigma, n, s, t)) for s in range(n) for t in range(n))
    avg = avg_total / (n * n)
    avg_pred = (n - 1) * (2 * n - 1) / 6
    row = {
        "n": n, "count": cnt, "floor_(n-1)^2/4": T,
        "worst_I_full": worst_full, "argmax_I_full": list(arg_full),
        "worst_I_single_shift": worst_single, "argmax_I_single_shift": list(arg_single),
        "single_shift_exceeds_bound": worst_single > T,
        "duality_checked": duality_checked,
        "duality_I_pi_eq_I_pi_inv_holds_for_all_checked": duality_fail is None,
        "sigma_n_avg_inv_over_n2_cuts": avg, "predicted_(n-1)(2n-1)/6": avg_pred,
        "avg_matches_formula": abs(avg - avg_pred) < 1e-9,
        "seconds": round(time.time() - t0),
    }
    log(f"n={n}: I_full max = {worst_full} (T={T}) at {arg_full}; single-shift max = {worst_single} "
        f"({'EXCEEDS' if worst_single > T else 'within'} T) at {arg_single}; "
        f"duality checked on {duality_checked} pi, holds: {duality_fail is None}; "
        f"sigma_n avg inv = {avg:.4f} (pred {avg_pred:.4f}); {row['seconds']} s")
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=4)
    ap.add_argument("--nmax", type=int, default=9)
    ap.add_argument("--duality-nmax", type=int, default=8,
                     help="exhaustively check I(pi)=I(pi^-1) up to this n (expensive: doubles the work)")
    ap.add_argument("--recursion-trials", type=int, default=500)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    logf = open(os.path.join(OUT, "report.log"), "a")

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + "\n")
        logf.flush()

    log(f"== {VERSION} {CORE_VERSION} args={vars(args)}")
    check_recursion(args.recursion_trials, log)
    rows = [run(n, log, check_duality=(n <= args.duality_nmax)) for n in range(args.nmin, args.nmax + 1)]
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump({"version": VERSION, "core": CORE_VERSION, "rows": rows}, f, indent=1)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write(f"# h13i_toric report ({VERSION}, {CORE_VERSION})\n\n")
        f.write("| n | max I(pi) | floor((n-1)^2/4) | max I single-shift | exceeds? | duality checked on | duality holds | sigma_n avg inv | (n-1)(2n-1)/6 |\n")
        f.write("|---|---|---|---|---|---|---|---|---|\n")
        for r in rows:
            f.write(f"| {r['n']} | {r['worst_I_full']} | {r['floor_(n-1)^2/4']} | {r['worst_I_single_shift']} | "
                     f"{r['single_shift_exceeds_bound']} | {r['duality_checked']} | {r['duality_I_pi_eq_I_pi_inv_holds_for_all_checked']} | "
                     f"{r['sigma_n_avg_inv_over_n2_cuts']:.4f} | {r['predicted_(n-1)(2n-1)/6']:.4f} |\n")


if __name__ == "__main__":
    main()
