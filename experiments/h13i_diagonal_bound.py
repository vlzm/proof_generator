"""H13-I attempt: does averaging pre-minimised "diagonal" cuts prove I(pi) <= floor((n-1)^2/4)?

Recap (docs/notes/h13_line_model.md): a double cut (q, c) gives the line
w_j = pi(q+1+j) - (q+1-c) mod n, j = 0..n-1, and I(pi) = min_{q,c} inv(w) over
all n^2 cuts. H13-I claims I(pi) <= floor((n-1)^2/4) for every pi, n >= 4.
Plain averaging of inv(w) over all n^2 cuts is already known NOT to prove this
(docs/notes/h13_line_model.md section 1).

New reduction tried in this session. Write k = (c - q) mod n. A short
computation (see write-up in h13_line_model.md section 7) shows that moving
along one "diagonal" of the (q, c) grid -- i.e. fixing k and letting q range
over 0..n-1 -- keeps the value-shift term (q+1-c) constant (it equals 1-k).
So on diagonal k, the n lines obtained are exactly the n cyclic ROTATIONS
(no value relabelling) of ONE fixed array

    u^(k)_i = (pi(i) + k - 1) mod n,   i = 0..n-1.

Hence, writing R(u) = min over the n cyclic rotations of u of inv(rotation),

    I(pi) = min_k R(u^(k))   (exact identity, cross-checked below against the
                               brute-force n^2-cut definition)

and in particular

    I(pi) <= (1/n) * sum_{k=0}^{n-1} R(u^(k))                              (*)

Candidate claim tested here:

    H13-I-diag:  sum_{k=0}^{n-1} R(u^(k)) <= n * floor((n-1)^2/4)   for all pi, n >= 4,

with equality iff pi is (a rotation of) a reflection pi_h(i) = h - i mod n.
If true this would prove H13-I via (*). It is STRICTLY STRONGER than H13-I
(average of the n diagonal minima vs. the true minimum), so it is a genuine,
falsifiable, more specific target -- not a restatement.

Status reached this session: VERIFIED, no counterexample, at
  - exhaustive enumeration of all pi, 4 <= n <= 10 (--exhaustive, this file);
  - structured/random families at n = 20, 50, 100, 101, 200 (--sample):
    sigma_n, reflection h=0, rev_n, id_n, affine maps, 3 random pi each.
Equality (sum == n * floor((n-1)^2/4) exactly) holds on every reflection
tested and nowhere else observed.

NOT proved: no closed-form argument for H13-I-diag was found this session.
A natural weaker route -- replacing R(u^(k)) (true minimum over rotations)
by inv(u^(k)) at a SPECIFIC explicit cut, e.g. "cut right after the position
holding the maximum value" -- fails already at n = 4 (see --heuristic mode,
and the negative result logged in the report): the true rotation-minimum is
essential, an explicit one-line rule is not enough. A second natural route --
bounding the minimum by the AVERAGE over rotations t of the same walk -- is
also too weak (see --raw mode): the un-rotated / averaged sum exceeds
n * floor((n-1)^2/4) already at n = 7 (91 vs 63 for sigma_7). Only the true
per-k minimum, then averaged over k, closes the gap in every case tested.

One auxiliary fact IS proved (elementary, not just checked): for the rotation
walk S_t^(k) = sum_{s<t} ((n-1) - 2 u^(k)_s) (cumulative inversion change when
rotating u^(k) by t), we have sum_{k=0}^{n-1} S_t^(k) = 0 for every fixed t.
Proof: for fixed i, as k ranges over 0..n-1, u^(k)_i = (pi(i)+k-1) mod n takes
every value in {0,...,n-1} exactly once, so sum_k u^(k)_i = n(n-1)/2 and
sum_k [(n-1) - 2 u^(k)_i] = n(n-1) - n(n-1) = 0; sum over i < t preserves this.
This identity is checked in --identity mode. It is NOT by itself enough to
prove H13-I-diag (see --raw mode above): using it via "min_t <= avg_t" loses
too much.

Usage:
  python3 experiments/h13i_diagonal_bound.py --exhaustive --nmin 4 --nmax 10
  python3 experiments/h13i_diagonal_bound.py --sample --sizes 20,50,100,101,200
  python3 experiments/h13i_diagonal_bound.py --heuristic --nmin 4 --nmax 8
  python3 experiments/h13i_diagonal_bound.py --raw --sizes 7,8,20
  python3 experiments/h13i_diagonal_bound.py --identity --nmin 4 --nmax 7

Version h13i_diagonal_bound-1.0 (session 9).
"""

import argparse
import itertools
import json
import os
import random
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def inv_count(arr):
    n = len(arr)
    c = 0
    for i in range(n):
        ai = arr[i]
        for j in range(i + 1, n):
            if ai > arr[j]:
                c += 1
    return c


def bound_I(n):
    """floor((n-1)^2/4), the conjectured max of I(pi) (C33/H13-I)."""
    return ((n - 1) ** 2) // 4


def line_inv(pi, q, c):
    """Direct definition (docs/notes/h13_line_model.md section 0): brute-force I(pi)."""
    n = len(pi)
    off = (q + 1 - c) % n
    w = [(pi[(q + 1 + j) % n] - off) % n for j in range(n)]
    return inv_count(w)


def I_pi_bruteforce(pi):
    n = len(pi)
    best = None
    for q in range(n):
        for c in range(n):
            iv = line_inv(pi, q, c)
            if best is None or iv < best:
                best = iv
    return best


def R_of_rotations(u):
    """min over the n cyclic rotations of u (values fixed, positions rotated) of inv.

    Incremental step (proved in the write-up): moving the front element with
    value m to the back changes inv by (n-1) - 2m.
    """
    n = len(u)
    cur = list(u)
    iv = inv_count(cur)
    best = iv
    for _ in range(1, n):
        m = cur[0]
        iv = iv + (n - 1) - 2 * m
        cur = cur[1:] + cur[:1]
        if iv < best:
            best = iv
    return best


def diagonal_values(pi):
    """R(u^(k)) for k = 0..n-1, u^(k)_i = (pi(i) + k - 1) mod n."""
    n = len(pi)
    return [R_of_rotations([(v + (k - 1) % n) % n for v in pi]) for k in range(n)]


def I_via_diagonals(pi):
    return min(diagonal_values(pi))


def sum_diagonals(pi):
    return sum(diagonal_values(pi))


def cut_after_max_heuristic(pi):
    """Explicit (non-optimal) per-diagonal cut: rotate to put the max value last.
    Tested as a candidate replacement for the true minimum; see module docstring.
    """
    n = len(pi)
    total = 0
    for k in range(n):
        u = [(v + (k - 1) % n) % n for v in pi]
        p = u.index(n - 1)
        t = (p + 1) % n
        arr = u[t:] + u[:t]
        total += inv_count(arr)
    return total


def raw_sum_no_rotation(pi):
    """sum_k inv(u^(k)) at t = 0 (no rotation at all); tests the weaker
    'average over t' route, which module docstring shows is insufficient.
    """
    n = len(pi)
    total = 0
    for k in range(n):
        u = [(v + (k - 1) % n) % n for v in pi]
        total += inv_count(u)
    return total


def walk_S(u):
    """S_t for t = 0..n (S_0 = S_n = 0), steps (n-1) - 2*u_s."""
    n = len(u)
    S = [0] * (n + 1)
    for t in range(1, n + 1):
        S[t] = S[t - 1] + (n - 1) - 2 * u[t - 1]
    return S


def cmd_exhaustive(nmin, nmax, verify_identity, out):
    report = {"mode": "exhaustive", "results": []}
    for n in range(nmin, nmax + 1):
        t0 = time.time()
        b = bound_I(n)
        worst_sum = -1
        worst_pi = None
        fails = 0
        cnt = 0
        equal_count = 0
        cross_check_n = min(n, 7)  # brute-force I(pi) cross-check only for small n (cheap)
        for pi in itertools.permutations(range(n)):
            cnt += 1
            pl = list(pi)
            s = sum_diagonals(pl)
            if s > n * b:
                fails += 1
            if s == n * b:
                equal_count += 1
            if s > worst_sum:
                worst_sum = s
                worst_pi = pi
            if n <= cross_check_n:
                assert I_via_diagonals(pl) == I_pi_bruteforce(pl), (
                    "diagonal decomposition disagrees with brute-force I(pi)",
                    pl,
                )
        dt = time.time() - t0
        line = {
            "n": n,
            "count": cnt,
            "bound": b,
            "n_times_bound": n * b,
            "fails": fails,
            "equal_to_bound_count": equal_count,
            "worst_sum": worst_sum,
            "worst_pi": list(worst_pi),
            "seconds": round(dt, 1),
        }
        report["results"].append(line)
        print(
            f"n={n}: count={cnt} bound={b} n*bound={n*b} fails={fails} "
            f"equal={equal_count} worst_sum={worst_sum} worst_pi={worst_pi} "
            f"time={dt:.1f}s"
        )
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as f:
        json.dump(report, f, indent=2)
    print("Saved", out)


def cmd_sample(sizes, out):
    report = {"mode": "sample", "results": []}
    rng = random.Random(42)
    for n in sizes:
        b = bound_I(n)
        families = {}
        families["sigma_n (h=1)"] = [(1 - i) % n for i in range(n)]
        families["reflection h=0"] = [(0 - i) % n for i in range(n)]
        families["rev_n"] = list(range(n - 1, -1, -1))
        families["id_n"] = list(range(n))
        a = 7
        while __import__("math").gcd(a, n) != 1:
            a += 1
        families[f"affine a={a},b=11"] = [(a * i + 11) % n for i in range(n)]
        for trial in range(3):
            p = list(range(n))
            rng.shuffle(p)
            families[f"random#{trial}"] = p
        for label, pi in families.items():
            s = sum_diagonals(pi)
            ok = s <= n * b
            report["results"].append(
                {"n": n, "label": label, "sum": s, "n_times_bound": n * b, "bound": b, "ok": ok}
            )
            print(
                f"n={n:4d} {label:20s} sum={s:9d} n*bound={n*b:9d} "
                f"avg={s/n:10.3f} bound={b:7d} {'OK' if ok else 'FAIL'}"
            )
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as f:
        json.dump(report, f, indent=2)
    print("Saved", out)


def cmd_heuristic(nmin, nmax, out):
    """Tests whether the explicit 'cut after max value' rule can replace the
    true rotation-minimum. Expected (and confirmed): it fails already at n=4.
    """
    report = {"mode": "heuristic_cut_after_max", "results": []}
    for n in range(nmin, nmax + 1):
        b = bound_I(n)
        fails = 0
        total = 0
        worst = -1
        worst_pi = None
        for pi in itertools.permutations(range(n)):
            total += 1
            s = cut_after_max_heuristic(list(pi))
            if s > n * b:
                fails += 1
            if s > worst:
                worst = s
                worst_pi = pi
        line = {"n": n, "bound": b, "n_times_bound": n * b, "fails": fails, "total": total, "worst": worst}
        report["results"].append(line)
        print(f"n={n}: bound={b} n*bound={n*b} fails={fails}/{total} worst={worst} at {worst_pi}")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as f:
        json.dump(report, f, indent=2)
    print("Saved", out)


def cmd_raw(sizes, out):
    """Tests the weaker 'average over t, not minimum' route. Expected (and
    confirmed): sum_k inv(u^(k)) at t=0 alone already exceeds n*bound on
    reflections for n >= 7, so this route cannot work either.
    """
    report = {"mode": "raw_no_rotation", "results": []}
    for n in sizes:
        b = bound_I(n)
        refl = [(1 - i) % n for i in range(n)]
        s = raw_sum_no_rotation(refl)
        report["results"].append({"n": n, "sum": s, "n_times_bound": n * b})
        print(f"n={n}: sum_k inv(u^k) at t=0 (sigma_n) = {s}, n*bound = {n*b}, diff={s-n*b}")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as f:
        json.dump(report, f, indent=2)
    print("Saved", out)


def cmd_identity(nmin, nmax, out):
    """Checks sum_k S_t^(k) == 0 for every fixed t, all pi, all t (proved fact)."""
    report = {"mode": "identity_sum_S_t", "results": []}
    for n in range(nmin, nmax + 1):
        bad = 0
        checked = 0
        for pi in itertools.permutations(range(n)):
            pl = list(pi)
            totals = [0] * (n + 1)
            for k in range(n):
                u = [(v + (k - 1) % n) % n for v in pl]
                S = walk_S(u)
                for t in range(n + 1):
                    totals[t] += S[t]
            checked += 1
            if any(x != 0 for x in totals):
                bad += 1
        report["results"].append({"n": n, "checked": checked, "bad": bad})
        print(f"n={n}: checked {checked} permutations, identity violated on {bad}")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as f:
        json.dump(report, f, indent=2)
    print("Saved", out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exhaustive", action="store_true")
    ap.add_argument("--sample", action="store_true")
    ap.add_argument("--heuristic", action="store_true")
    ap.add_argument("--raw", action="store_true")
    ap.add_argument("--identity", action="store_true")
    ap.add_argument("--nmin", type=int, default=4)
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--sizes", type=str, default="20,50,100,101,200")
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    out_dir = os.path.join(ROOT, "data", "runs", "h13i_diagonal_bound")
    if args.exhaustive:
        out = args.out or os.path.join(out_dir, f"exhaustive_{args.nmin}_{args.nmax}.json")
        cmd_exhaustive(args.nmin, args.nmax, True, out)
    elif args.sample:
        sizes = [int(x) for x in args.sizes.split(",")]
        out = args.out or os.path.join(out_dir, "sample.json")
        cmd_sample(sizes, out)
    elif args.heuristic:
        out = args.out or os.path.join(out_dir, f"heuristic_{args.nmin}_{args.nmax}.json")
        cmd_heuristic(args.nmin, args.nmax, out)
    elif args.raw:
        sizes = [int(x) for x in args.sizes.split(",")]
        out = args.out or os.path.join(out_dir, "raw.json")
        cmd_raw(sizes, out)
    elif args.identity:
        out = args.out or os.path.join(out_dir, f"identity_{args.nmin}_{args.nmax}.json")
        cmd_identity(args.nmin, args.nmax, out)
    else:
        print("Choose one of --exhaustive / --sample / --heuristic / --raw / --identity")
        sys.exit(1)


if __name__ == "__main__":
    main()
