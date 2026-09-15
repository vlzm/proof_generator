"""H13-I proof attempt (session 9): rotation-recurrence lemma + candidate cuts.

H13-I: for every permutation pi of Z_n (n >= 4) there is a double cut (a, b)
(position shift a, value shift b) with inv(w^{a,b}) <= floor((n-1)^2/4), where
w^{a,b}_j = (pi((a+j) mod n) - b) mod n, j = 0..n-1.  Equivalently (rule: see
docs/notes/h13_line_model.md #7) max_{a,b} NonInv(a,b) >= floor(n^2/4), where
NonInv(a,b) = C(n,2) - inv(w^{a,b}).

This script checks (independently of oracle/moves.py -- direct inversion
counting only, per the task) the two pieces of this session's investigation:

1. Lemma D (rotation recurrence) and its dual: for fixed b,
     NonInv(a+1,b) - NonInv(a,b) = 2*r_a - (n-1),   r_a = (pi(a)-b) mod n;
   for fixed a,
     NonInv(a,b+1) - NonInv(a,b) = 2*r'_b - (n-1),  r'_b = (pi^{-1}(b)-a) mod n.
   Both are proved directly (docs/notes/h13_line_model.md #7, Lemma D); this
   verifies them by brute force, exhaustively for 4 <= n <= NMAX and by sample
   for larger n.

2. Five candidate O(n)-size families of cuts, each tested for whether
   min over the family of inv(w^{a,b}) already meets floor((n-1)^2/4) for
   every pi (exhaustive n = 4..NMAX):
     - diag:      (a, pi(a)) for a in Z_n   (w_0 = 0 always)
     - a0:        (0, b) for b in Z_n       (position cut fixed at 0)
     - b0:        (a, 0) for a in Z_n       (value cut fixed at 0)
     - a0_or_b0:  union of the two above
     - abdiag:    (t, t) for t in Z_n
   Also the second-moment (RMS) bound: is sqrt(mean_{a,b} NonInv(a,b)^2)
   already >= floor(n^2/4) for every pi?  (A positive answer would prove
   H13-I via Cauchy-Schwarz max >= RMS; checked exhaustively n = 4..NMAX.)

All are expected (and found) to fail for some pi at small n -- this is a
negative/exploratory experiment, not a proof.  See the note for the verdict.

Usage: python3 experiments/h13i_walk_attempt.py --nmax 8
"""

import argparse
import itertools
import json
import math
import random
import time


def inv_count(seq):
    n = len(seq)
    c = 0
    for j in range(n):
        sj = seq[j]
        for k in range(j + 1, n):
            if sj > seq[k]:
                c += 1
    return c


def noninv(pi, a, b, n, C):
    w = [(pi[(a + j) % n] - b) % n for j in range(n)]
    return C - inv_count(w)


def all_noninv(pi, n, C):
    return [noninv(pi, a, b, n, C) for a in range(n) for b in range(n)]


def check_recurrence(pi, n):
    """Verify both directions of Lemma D exactly. Returns True/False."""
    inv_pi = [0] * n
    for i, v in enumerate(pi):
        inv_pi[v] = i
    C = n * (n - 1) // 2
    for b in range(n):
        vals = [noninv(pi, a, b, n, C) for a in range(n)]
        for a in range(n):
            lhs = vals[(a + 1) % n] - vals[a]
            r_a = (pi[a] - b) % n
            if lhs != 2 * r_a - (n - 1):
                return False, ("a", a, b)
    for a in range(n):
        vals = [noninv(pi, a, b, n, C) for b in range(n)]
        for b in range(n):
            lhs = vals[(b + 1) % n] - vals[b]
            rp_b = (inv_pi[b] - a) % n
            if lhs != 2 * rp_b - (n - 1):
                return False, ("b", a, b)
    return True, None


def family_min(pi, n, kind, C):
    """min_{(a,b) in family} inv(w^{a,b}) = C - max_{(a,b) in family} NonInv(a,b)."""
    if kind == "diag":
        return C - max(noninv(pi, a, pi[a], n, C) for a in range(n))
    if kind == "a0":
        return C - max(noninv(pi, 0, b, n, C) for b in range(n))
    if kind == "b0":
        return C - max(noninv(pi, a, 0, n, C) for a in range(n))
    if kind == "a0_or_b0":
        return C - max(
            max(noninv(pi, 0, b, n, C) for b in range(n)),
            max(noninv(pi, a, 0, n, C) for a in range(n)),
        )
    if kind == "abdiag":
        return C - max(noninv(pi, t, t, n, C) for t in range(n))
    raise ValueError(kind)


FAMILIES = ["diag", "a0", "b0", "a0_or_b0", "abdiag"]


def run(nmax, nmin=4, seed=0, sample_large=(9, 10, 12, 15)):
    report = {"recurrence": {}, "families": {}, "rms": {}, "sample_large_recurrence": {}}

    rng = random.Random(seed)

    for n in range(nmin, nmax + 1):
        t0 = time.time()
        C = n * (n - 1) // 2
        bound_I = (n - 1) ** 2 // 4
        bound_R = n * n // 4

        rec_fail = None
        family_worst = {k: (-1, None) for k in FAMILIES}
        rms_worst = (1e18, None)

        count = 0
        for pi in itertools.permutations(range(n)):
            count += 1
            ok, where = check_recurrence(pi, n)
            if not ok:
                rec_fail = (pi, where)

            for k in FAMILIES:
                m = family_min(pi, n, k, C)
                if m > family_worst[k][0]:
                    family_worst[k] = (m, pi)

            vals = all_noninv(pi, n, C)
            mean2 = sum(v * v for v in vals) / len(vals)
            rms = math.sqrt(mean2)
            margin = rms - bound_R
            if margin < rms_worst[0]:
                rms_worst = (margin, pi)

        report["recurrence"][n] = "OK" if rec_fail is None else f"FAIL {rec_fail}"
        report["families"][n] = {
            "bound_I": bound_I,
            "per_family_max_min": {k: family_worst[k][0] for k in FAMILIES},
            "per_family_worst_pi": {k: family_worst[k][1] for k in FAMILIES},
            "all_families_sufficient": all(family_worst[k][0] <= bound_I for k in FAMILIES),
        }
        report["rms"][n] = {
            "bound_R": bound_R,
            "worst_margin": rms_worst[0],
            "worst_pi": rms_worst[1],
            "rms_bound_sufficient": rms_worst[0] >= 0,
        }
        print(
            f"n={n} perms={count} time={time.time()-t0:.1f}s "
            f"recurrence={report['recurrence'][n]} "
            f"families_ok={report['families'][n]['all_families_sufficient']} "
            f"rms_ok={report['rms'][n]['rms_bound_sufficient']}"
        )

    for n in sample_large:
        fails = 0
        trials = 30
        for _ in range(trials):
            pi = list(range(n))
            rng.shuffle(pi)
            ok, where = check_recurrence(tuple(pi), n)
            if not ok:
                fails += 1
        report["sample_large_recurrence"][n] = f"{trials - fails}/{trials} OK"
        print(f"n={n} (sampled {trials} random perms) recurrence: {trials-fails}/{trials} OK")

    return report


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=4)
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=str, default="data/runs/h13i_attempt/report.json")
    args = ap.parse_args()

    rep = run(args.nmax, args.nmin, args.seed)

    import os

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(rep, f, indent=2, default=str)
    print("wrote", args.out)
