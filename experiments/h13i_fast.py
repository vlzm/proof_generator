"""H13-I: exact O(n^2) formula for I(pi) and extended numerical verification
(session 9).

Derivation (see docs/proofs/C37_value_shift_formula.md for the write-up and
verification). For a position cut a and value cut b (line model notation,
docs/notes/h13_line_model.md §0; a = q+1, b = q+1-c), inv(a,b) as a function
of b (a fixed) satisfies the exact telescoping identity

    inv(a, b+1) - inv(a, b) = (n-1) - 2*pos_{a,b}(0),

where pos_{a,b}(0) is the domain index j with pi((a+j) mod n) = b (the point
currently mapped to the value being "peeled off"). Writing sigma_a(t) =
(pi^{-1}(t) - a) mod n (so sigma_a is v^{(a)-1}, the inverse of the
a-rotation of pi) and inv(a, 0) = inv_count(v^{(a)}) = inv_count(sigma_a)
(inversions of a permutation equal inversions of its inverse), this gives a
closed form for the whole b-optimization:

    min_b inv(a, b) = inv_count(sigma_a) - max_{0<=b<=n} [2*S_a(b) - b(n-1)],
    S_a(b) = sigma_a(0) + ... + sigma_a(b-1).

So I(pi) = min_a of the right-hand side, an O(n) evaluation per a (O(n log n)
for inv_count once, then inv_count(sigma_{a+1}) = inv_count(sigma_a) +
(n-1) - 2*pi(a), an O(1) step) -- O(n^2) total per pi instead of the O(n^4)
brute n^2-cut x n^2-inversion-count search. This does not by itself prove
H13-I (the b-side alone is not enough, docs/notes/h13_line_model.md Lemma C /
this session's G-universal test), but it lets the numerical ladder (AGENTS.md
rule 9) go much further than the O(n^4) brute force used for C33/C35.

Verified against the brute O(n^4) search on 200+ random permutations,
4 <= n <= 9 (exact match every time).

Usage:
  python3 experiments/h13i_fast.py --exhaustive 11          # extends C33 range
  python3 experiments/h13i_fast.py --sample 500 --n 100,500,1000,5000 --seed 1
Output: data/runs/h13i_fast/. Version h13i_fast-1.0.
"""

import argparse
import itertools
import json
import os
import random
import time

VERSION = "h13i_fast-1.0"
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "data", "runs", "h13i_fast")


def inv_count_bit(seq):
    n = len(seq)
    bit = [0] * (n + 1)

    def upd(i):
        i += 1
        while i <= n:
            bit[i] += 1
            i += i & (-i)

    def qry(i):
        i += 1
        s = 0
        while i > 0:
            s += bit[i]
            i -= i & (-i)
        return s

    inv = 0
    for x in reversed(seq):
        if x > 0:
            inv += qry(x - 1)
        upd(x)
    return inv


def fast_I(pi, n):
    """I(pi) = min_{a,b} inv(a,b), exact, O(n^2)."""
    pinv = [0] * n
    for i, v in enumerate(pi):
        pinv[v] = i
    sigma0 = pinv
    inv_a = inv_count_bit(sigma0)
    invs = [0] * n
    for a in range(n):
        invs[a] = inv_a
        inv_a = inv_a + (n - 1) - 2 * pi[a]
    best = None
    argbest = None
    for a in range(n):
        S = 0
        bestE = 0
        for t in range(n):
            S += 2 * ((sigma0[t] - a) % n) - (n - 1)
            if S > bestE:
                bestE = S
        Ga = invs[a] - bestE
        if best is None or Ga < best:
            best, argbest = Ga, a
    return best, argbest


def run_exhaustive(n, checkpoint_every=2_000_000):
    bound = ((n - 1) ** 2) // 4
    t0 = time.time()
    worst = -1
    worst_pi = None
    count = 0
    for pi in itertools.permutations(range(n)):
        I, _ = fast_I(list(pi), n)
        if I > worst:
            worst, worst_pi = I, pi
        count += 1
        if count % checkpoint_every == 0:
            print(f"  ... {count} perms, elapsed {time.time()-t0:.1f}s, "
                  f"worst so far {worst}")
    dt = time.time() - t0
    return {
        "mode": "exhaustive", "n": n, "bound": bound, "count": count,
        "time_s": round(dt, 1), "max_I": worst, "max_I_witness": list(worst_pi),
        "ok": worst <= bound,
    }


def run_sample(n, sample, seed):
    bound = ((n - 1) ** 2) // 4
    random.seed(seed)
    base = list(range(n))
    t0 = time.time()
    worst = -1
    worst_pi = None
    for _ in range(sample):
        random.shuffle(base)
        I, _ = fast_I(base, n)
        if I > worst:
            worst, worst_pi = I, list(base)
    dt = time.time() - t0
    return {
        "mode": "sampled", "n": n, "bound": bound, "count": sample, "seed": seed,
        "time_s": round(dt, 2), "max_I": worst, "max_I_witness": worst_pi,
        "ok": worst <= bound,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exhaustive", type=int, default=None,
                     help="run exhaustive search up to and including this n "
                          "(starting from 4)")
    ap.add_argument("--sample", type=int, default=0)
    ap.add_argument("--n", type=str, default="")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    results = {"version": VERSION, "started": time.strftime("%Y-%m-%d %H:%M:%S"),
               "runs": []}

    if args.exhaustive:
        for n in range(4, args.exhaustive + 1):
            r = run_exhaustive(n)
            results["runs"].append(r)
            print(f"[exhaustive] n={n} count={r['count']} time={r['time_s']}s "
                  f"max_I={r['max_I']} bound={r['bound']} ok={r['ok']}")

    if args.sample and args.n:
        for n in [int(x) for x in args.n.split(",")]:
            r = run_sample(n, args.sample, args.seed)
            results["runs"].append(r)
            print(f"[sampled] n={n} count={r['count']} time={r['time_s']}s "
                  f"max_I={r['max_I']} bound={r['bound']} ok={r['ok']}")

    with open(os.path.join(OUT, "report.json"), "w") as fh:
        json.dump(results, fh, indent=2)
    print("Report written to", OUT)


if __name__ == "__main__":
    main()
