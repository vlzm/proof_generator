"""Arc decomposition of the double-cut inversion count (H13-I, session 9).

For pi in S_n and a double cut (q, c) (PLAN/h13_line_model.md conventions:
w_j = pi(q+1+j) - (q+1-c) mod n), write t = (q+1-c) mod n (the value shift;
(q, c) <-> (q, t) is a bijection of Z_n x Z_n for fixed q). For each pair of
positions i < i' (as naturals in [0, n)) let L = i'-i, C = 1 if pi(i) > pi(i')
else 0, vmin = min(pi(i), pi(i')), M = |pi(i)-pi(i')|. Define the two
arc-indicators (both independent of pi's *values* resp. *positions*):
  A(q)  = 1 iff q in [i, i+L-1]           (position order of i, i' reverses)
  R(t)  = 1 iff t in [vmin+1, vmin+M]     (value order of the pair reverses)
Then the pair is an inversion of the line iff A(q) xor C xor R(t) = 1, so
    inv(q, t) = N1(q) + beta(t) - 2*Gamma(q, t)
where
    N1(q)      = #{pairs : A(q) xor C = 1}            (depends on pi and q)
    beta(t)    = t*(n-t)                              (universal, independent
                                                          of pi -- see below)
    Gamma(q,t) = #{pairs with A(q) xor C = 1 and R(t) = 1}.

beta(t) is universal because the multiset of unordered value-pairs {pi(i),
pi(i')} over all position-pairs is exactly the set of all C(n,2) pairs of
Z_n (pi is a bijection), so summing R(t) over all pairs counts, for each t,
the number of value-pairs whose vmin < t <= vmax, i.e. t*(n-t). (Similarly
alpha(q) = sum_pairs A(q) = (q+1)*(n-1-q) is universal, but does not appear
in the decomposition above since N1(q) mixes A(q) with the pi-dependent bit
C via xor.)

This gives an *exact* identity (checked against brute force for all pi at
4 <= n <= 6, exhaustively) refining Lemma C of docs/notes/h13_line_model.md
(which only gave linear *inequalities* satisfied at a local optimum) into a
closed form for every cut. Two consequences, both checked here:

1. Gamma(q,t) <= min(N1(q), beta(t)) trivially, hence
   inv(q, t) >= |N1(q) - beta(t)|  for every (q, t).
   This is a *necessary* condition on any (q, t) that could achieve
   I(pi) <= floor((n-1)^2/4): N1(q) and beta(t) must be close.

2. The natural "coupled shift" candidate from h13_line_model.md section 6 --
   choose q minimizing N1(q) alone (ignoring the value side), then choose
   the best t for that q -- is tested exhaustively and FAILS already at
   n = 7 (witness recorded below): it does not by itself yield
   inv <= floor((n-1)^2/4) for every pi, even though the true minimum
   I(pi) = min_{q,t} inv(q,t) does satisfy the bound at this n (C33/C35).
   So minimizing the two "sides" of the decomposition independently is not
   enough; the interaction term Gamma(q,t) must be controlled jointly, as
   already observed for plain averaging (h13_line_model.md section 1.7).

Usage: python3 experiments/h13i_arc_decomposition.py --nmax 8
Output: data/runs/h13i_arc_decomposition/. Version h13i_arc_decomposition-1.0.
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

VERSION = "h13i_arc_decomposition-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h13i_arc_decomposition")


def floor_half_sq(n):
    return ((n - 1) ** 2) // 4


def inv_count(seq):
    n = len(seq)
    return sum(1 for i in range(n) for j in range(i + 1, n) if seq[i] > seq[j])


def line_for(pi, q, c, n):
    shift = (q + 1 - c) % n
    return [(pi[(q + 1 + j) % n] - shift) % n for j in range(n)]


def brute_I(pi, n):
    best = None
    for q in range(n):
        for c in range(n):
            iv = inv_count(line_for(pi, q, c, n))
            if best is None or iv < best:
                best = iv
    return best


def pair_data(pi, n):
    """(i, i', L, vmin, M, C) for every pair i < i' (naturals)."""
    out = []
    for i in range(n):
        for ip in range(i + 1, n):
            L = ip - i
            v, vp = pi[i], pi[ip]
            vmin, vmax = (v, vp) if v < vp else (vp, v)
            out.append((i, ip, L, vmin, vmax - vmin, 1 if v > vp else 0))
    return out


def A(q, i, L):
    return 1 if i <= q <= i + L - 1 else 0


def R(t, vmin, M):
    return 1 if vmin + 1 <= t <= vmin + M else 0


def beta(t, n):
    return t * (n - t)


def inv_via_decomposition(data, q, t):
    n1 = 0
    gamma = 0
    for (i, ip, L, vmin, M, C) in data:
        d = A(q, i, L) ^ C
        if d:
            n1 += 1
            if R(t, vmin, M):
                gamma += 1
    return n1, gamma


def check_identity(nmax, log):
    for n in range(4, nmax + 1):
        ok = True
        for pi in itertools.permutations(range(n)):
            data = pair_data(pi, n)
            for q in range(n):
                for c in range(n):
                    t = (q + 1 - c) % n
                    iv_direct = inv_count(line_for(pi, q, c, n))
                    n1, gamma = inv_via_decomposition(data, q, t)
                    iv_formula = n1 + beta(t, n) - 2 * gamma
                    if iv_direct != iv_formula:
                        log(f"MISMATCH n={n} pi={pi} q={q} c={c} t={t} "
                            f"direct={iv_direct} formula={iv_formula}")
                        ok = False
        log(f"n={n}: identity inv(q,t)=N1(q)+beta(t)-2*Gamma(q,t) "
            f"{'OK (exhaustive)' if ok else 'FAILED'}")


def strategy_min_N1_then_t(data, n):
    """Choose q minimizing N1(q) alone (no look-ahead on t), then the best t
    for that q. Tests the 'coupled shift' idea of h13_line_model.md section 6
    reduced to its simplest form: optimize the two sides of the decomposition
    separately instead of jointly."""
    n1_by_q = []
    for q in range(n):
        n1 = sum(1 for (i, ip, L, vmin, M, C) in data if (A(q, i, L) ^ C))
        n1_by_q.append(n1)
    qbest = min(range(n), key=lambda q: n1_by_q[q])
    n1min = n1_by_q[qbest]
    best_inv = None
    for t in range(n):
        gamma = sum(1 for (i, ip, L, vmin, M, C) in data
                    if (A(qbest, i, L) ^ C) and R(t, vmin, M))
        iv = n1min + beta(t, n) - 2 * gamma
        if best_inv is None or iv < best_inv:
            best_inv = iv
    return best_inv


def run_strategy_scan(nmax, log):
    rows = []
    for n in range(4, nmax + 1):
        t0 = time.time()
        thresh = floor_half_sq(n)
        max_I = 0
        max_strat = 0
        worst_strat = None
        fails = []
        for pi in itertools.permutations(range(n)):
            data = pair_data(pi, n)
            I = brute_I(pi, n)
            strat = strategy_min_N1_then_t(data, n)
            assert strat >= I, (pi, strat, I)
            if I > max_I:
                max_I = I
            if strat > max_strat:
                max_strat = strat
                worst_strat = pi
            if strat > thresh:
                fails.append(list(pi))
        elapsed = time.time() - t0
        log(f"n={n}: floor((n-1)^2/4)={thresh} max_I={max_I} "
            f"max(strategy min-N1-then-t)={max_strat} "
            f"worst_witness={worst_strat} fails={len(fails)} ({elapsed:.1f}s)")
        rows.append({
            "n": n, "threshold": thresh, "max_I": max_I,
            "max_strategy_min_N1_then_t": max_strat,
            "strategy_worst_witness": list(worst_strat) if worst_strat else None,
            "strategy_fail_count": len(fails),
            "strategy_fail_examples": fails[:5],
            "seconds": round(elapsed, 1),
        })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--identity_nmax", type=int, default=6,
                     help="exhaustive check of the exact identity, all (pi,q,c)")
    ap.add_argument("--nmax", type=int, default=8,
                     help="exhaustive strategy-vs-I comparison, all pi")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    logf = open(os.path.join(OUT, "report.log"), "a")

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + "\n")
        logf.flush()

    log(f"== {VERSION} {CORE_VERSION} args={vars(args)}")
    check_identity(args.identity_nmax, log)
    rows = run_strategy_scan(args.nmax, log)
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump({"version": VERSION, "core": CORE_VERSION, "rows": rows}, f, indent=1)


if __name__ == "__main__":
    main()
