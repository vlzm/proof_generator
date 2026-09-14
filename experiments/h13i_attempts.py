"""H13-I attempt (session 9): test two candidate reductions of the double-cut
search (q, c) in Z_n x Z_n to a single free parameter or to a small family,
against I(pi) <= floor((n-1)^2/4) (PLAN Sec.8, docs/notes/h13_line_model.md Sec.6).

I(pi) = min over (q, c) of inv(w), w_j = (pi[(q+1+j) mod n] - (q+1-c)) mod n,
j = 0..n-1 (same convention as experiments/line_model.py).

Candidate A ("single free parameter"): fix q = n-1 (i.e. read pi in its own
position order, no position rotation) and search only over c. Tests whether
position rotation is actually needed, or whether the value shift alone
already suffices for every pi (as a *fixed* linear arrangement of values).

Candidate B ("point-anchored cuts"): restrict the search to the n cuts with
w_0 = 0, i.e. for each i in 0..n-1, the cut that places the point (i, pi(i))
at the origin of the line (q = i-1, c = i - pi(i) mod n). n cuts instead of
n^2.

Both are exhaustive checks against full_I (all n^2 cuts) for 4 <= n <= 9.
Usage: python3 experiments/h13i_attempts.py --nmax 9
Version: h13i_attempts-1.0.
"""

import argparse
import itertools
import json
import os
import time

VERSION = "h13i_attempts-1.0"
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "data", "runs", "h13i_attempts")


def inv(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def line_at(pi, n, q, c):
    return tuple((pi[(q + 1 + j) % n] - (q + 1 - c)) % n for j in range(n))


def full_I(pi, n):
    best = None
    for q in range(n):
        for c in range(n):
            iv = inv(line_at(pi, n, q, c))
            if best is None or iv < best:
                best = iv
    return best


def single_param_I(pi, n):
    """Candidate A: q fixed (=n-1, i.e. natural position order), c free."""
    q = n - 1
    best = None
    for c in range(n):
        iv = inv(line_at(pi, n, q, c))
        if best is None or iv < best:
            best = iv
    return best


def anchored_I(pi, n):
    """Candidate B: only cuts with w_0 = 0 (n cuts, one per point at origin)."""
    best = None
    for i in range(n):
        q = (i - 1) % n
        c = (i - pi[i]) % n
        iv = inv(line_at(pi, n, q, c))
        if best is None or iv < best:
            best = iv
    return best


def run(n, log):
    t0 = time.time()
    B = (n - 1) ** 2 // 4
    worst_full = worst_A = worst_B = -1
    arg_full = arg_A = arg_B = None
    cnt = 0
    for pi in itertools.permutations(range(n)):
        cnt += 1
        fi = full_I(pi, n)
        ai = single_param_I(pi, n)
        bi = anchored_I(pi, n)
        if fi > worst_full:
            worst_full, arg_full = fi, pi
        if ai > worst_A:
            worst_A, arg_A = ai, pi
        if bi > worst_B:
            worst_B, arg_B = bi, pi
    row = {
        "n": n, "count": cnt, "floor_(n-1)^2/4": B, "seconds": round(time.time() - t0, 1),
        "full_I_max": worst_full, "full_I_argmax": list(arg_full),
        "single_param_I_max": worst_A, "single_param_I_argmax": list(arg_A),
        "anchored_I_max": worst_B, "anchored_I_argmax": list(arg_B),
    }
    log(f"n={n}: full max I={worst_full} (bound {B}, {'OK' if worst_full <= B else 'FAIL'}); "
        f"single-param (q fixed, c free) max={worst_A} at {arg_A} "
        f"({'OK' if worst_A <= B else 'FAIL, refutes candidate A'}); "
        f"anchored (n cuts, w_0=0) max={worst_B} at {arg_B} "
        f"({'OK' if worst_B <= B else 'FAIL, refutes candidate B'}); {row['seconds']} s")
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=4)
    ap.add_argument("--nmax", type=int, default=8)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    logf = open(os.path.join(OUT, "report.log"), "a")

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + "\n")
        logf.flush()

    log(f"== {VERSION} args={vars(args)}")
    rows = [run(n, log) for n in range(args.nmin, args.nmax + 1)]
    with open(os.path.join(OUT, f"report_n{args.nmin}_{args.nmax}.json"), "w") as f:
        json.dump({"version": VERSION, "rows": rows}, f, indent=1)


if __name__ == "__main__":
    main()
