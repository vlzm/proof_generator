"""H13-I attempt (session 9): reduction of I(pi) to a nested pair of closed walks,
and tests of several reduced-parameter families that were hoped to avoid the full
n^2 search over cuts (q, c).

Setup (docs/notes/h13_line_model.md #6): I(pi) = min_{a,b in Z_n} inv(w^{a,b}),
w^{a,b}_j = (pi((a+j) mod n) - b) mod n for j = 0..n-1 (a = q+1 is the position
cut, b = c the value shift, in the notation of experiments/line_profile.c).

Finding 1 (identity, PROVED, exact -- not just numerically checked): let
tau = pi^{-1}. Define, for a in Z_n, q_a(v) = (tau(v) - a) mod n (v = 0..n-1);
this is the line of *positions*, value-shifted by a -- i.e. a plays for tau
exactly the role b plays for pi. Then

    inv(q_a) = inv(pi) + a*(n-1) - 2*Q(a),    Q(a) = sum_{i<a} pi(i)          (W-A)

and, for the *original* two-parameter quantity, with P_a(k) = sum_{v<k} q_a(v):

    inv(w^{a,k}) = inv(q_a) + k*(n-1) - 2*P_a(k)                              (W-B)

Both are proved by a one-line telescoping argument: replacing a shift s by s+1
moves the element currently holding value 0 (in the s-shifted sequence) from
rank 0 to rank n-1 and decrements every other value by 1, which changes inv by
exactly (n-1) - 2*p, p = position of that element in the *fixed* index order.
(W-A) is exactly this argument applied to the sequence pi (shift variable a,
walking value 0,1,2,... away), instantiated at s=a-1; (W-B) is the same
argument applied to q_a (shift variable k). The pi <-> tau, a <-> b symmetry is
exact: (W-A) for pi's own values is literally the b-walk of pi read backwards.
Verified by direct computation against a naive O(n^2) inversion counter for
random permutations, 4 <= n <= 20 (function `verify_identities`).

Finding 2 (negative): every one-parameter (O(n)-size, instead of the full
n^2-size) family of cuts tried gives a *counterexample* pi with
min over the family of inv(w) > floor((n-1)^2/4):

  - fix a = 0, vary b only ("Lemma L"): FAILS at n = 7,
    pi = (0, 5, 4, 3, 2, 1, 6) (this is exactly the Lemma-C counterexample of
    docs/notes/h13_line_model.md #4, re-derived independently here as a
    byproduct of walk (W-B) with a = 0 fixed).
  - fix b = 0, vary a only (equivalent to the same statement for tau, by the
    pi <-> tau duality of Finding 1): FAILS, same shape of counterexample.
  - the diagonal families b = a and b = -a (n points each): FAIL already at
    n = 4.
  - the *union* of "a = 0, any b" and "b = 0, any a" (2n-1 points instead of
    n^2): still FAILS at n = 7, same permutation as Lemma L.

Conclusion for this session: I(pi) genuinely needs joint optimisation over
both parameters; no O(n)-size reduced family of cuts suffices. See
docs/notes/h13_line_model.md #7 for the write-up and the recommended next
step. This script only records the (exact) identity and the (exhaustive,
small n) negative results; it does not attempt the full n^2 x n! search
(that is experiments/line_profile.c / checks/check_C34.py part C, already run
for 4 <= n <= 10).

Usage: python3 experiments/h13_toric_reduction.py [--nmax 7] [--identity-nmax 20]
Output: data/runs/h13_toric_reduction/report.json, report.md.
Version h13_toric_reduction-1.0.
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

VERSION = "h13_toric_reduction-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h13_toric_reduction")


def inv(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def verify_identities(nmax, trials, log):
    """Exact check of (W-A) and (W-B) against a naive inversion counter."""
    rng = random.Random(20260916)
    for n in range(4, nmax + 1):
        for _ in range(trials):
            pi = list(range(n))
            rng.shuffle(pi)
            pi = tuple(pi)
            tau = [0] * n
            for i, v in enumerate(pi):
                tau[v] = i
            inv_pi = inv(pi)
            Q = [0] * (n + 1)
            for a in range(n):
                Q[a + 1] = Q[a] + pi[a]
            for a in range(n):
                q_a = [(tau[v] - a) % n for v in range(n)]
                inv_qa = inv(tuple(q_a))
                rhs_a = inv_pi + a * (n - 1) - 2 * Q[a]
                assert inv_qa == rhs_a, ("W-A", n, pi, a, inv_qa, rhs_a)
                Pa = [0] * (n + 1)
                for k in range(n):
                    Pa[k + 1] = Pa[k] + q_a[k]
                for k in range(n + 1):
                    w = tuple((pi[(a + j) % n] - k) % n for j in range(n))
                    lhs_b = inv(w)
                    rhs_b = inv_qa + k * (n - 1) - 2 * Pa[k]
                    assert lhs_b == rhs_b, ("W-B", n, pi, a, k, lhs_b, rhs_b)
        log(f"identities n={n}: (W-A) and (W-B) hold exactly on {trials} random perms, all a, all k")
    return {"nmax": nmax, "trials": trials, "status": "both identities verified exactly"}


def family_worst(n, cuts_fn):
    """cuts_fn(pi) -> iterable of (a, b); returns (worst_min, worst_pi)."""
    worst, worst_pi = -1, None
    for pi in itertools.permutations(range(n)):
        best = min(inv(tuple((pi[(a + j) % n] - b) % n for j in range(n))) for a, b in cuts_fn(pi, n))
        if best > worst:
            worst, worst_pi = best, pi
    return worst, worst_pi


def run_families(nmax, log):
    families = {
        "a=0, any b": lambda pi, n: [(0, b) for b in range(n)],
        "b=0, any a": lambda pi, n: [(a, 0) for a in range(n)],
        "b=a": lambda pi, n: [(a, a) for a in range(n)],
        "b=-a": lambda pi, n: [(a, (-a) % n) for a in range(n)],
        "union(a=0,any b ; b=0,any a)": lambda pi, n: [(0, b) for b in range(n)] + [(a, 0) for a in range(n)],
    }
    rows = []
    for name, fn in families.items():
        first_fail = None
        for n in range(4, nmax + 1):
            bound = (n - 1) ** 2 // 4
            t0 = time.time()
            worst, ex = family_worst(n, fn)
            ok = worst <= bound
            log(f"family [{name}] n={n}: worst(min over family) = {worst}, "
                f"floor((n-1)^2/4) = {bound} -> {'OK' if ok else 'FAIL'}"
                f"{'' if ok else f', counterexample {ex}'} ({time.time() - t0:.1f} s)")
            rows.append({"family": name, "n": n, "worst_min": worst, "bound": bound, "ok": ok,
                         "counterexample": list(ex) if not ok and ex is not None else None})
            if not ok and first_fail is None:
                first_fail = n
        log(f"family [{name}]: {'never fails up to n=' + str(nmax) if first_fail is None else 'first fails at n=' + str(first_fail)}")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=7, help="max n for the exhaustive family tests (7 is enough to see all failures)")
    ap.add_argument("--identity-nmax", type=int, default=20)
    ap.add_argument("--identity-trials", type=int, default=80)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    lines = []

    def log(msg):
        print(msg, flush=True)
        lines.append(msg)

    log(f"== {VERSION} {CORE_VERSION} args={vars(args)}")
    ident = verify_identities(args.identity_nmax, args.identity_trials, log)
    rows = run_families(args.nmax, log)

    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump({"version": VERSION, "core": CORE_VERSION, "args": vars(args),
                   "identities": ident, "families": rows}, f, indent=1)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# h13_toric_reduction — нестед-волк тождество и провал редуцированных семейств разрезов\n\n")
        f.write(f"Команда: `python3 experiments/h13_toric_reduction.py --nmax {args.nmax} "
                f"--identity-nmax {args.identity_nmax}`. Версия: {VERSION}, {CORE_VERSION}.\n\n```text\n")
        f.write("\n".join(lines) + "\n```\n")
    log("done")


if __name__ == "__main__":
    main()
