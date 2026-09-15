"""H13-I exploration (session 9): does rotating ONE axis alone suffice?

H13-I claims every permutation pi has a "double cut" (independent rotation of
positions by a and values by b) with <= floor((n-1)^2/4) inversions, i.e. its
torus orbit under Z_n x Z_n (row rotation x column rotation) contains a member
with >= floor(n^2/4) non-inversions ("concordant pairs").

This script tests the natural weaker sub-claim: fix the position order (a = 0,
i.e. use pi's own sequence unshifted) and only search over the value rotation
b.  If this 1-parameter search already reached floor(n^2/4) concordant pairs
for every pi, H13-I would follow trivially without ever needing to rotate
positions.  Exhaustive result: it does NOT (first failure at n = 7), so the
double cut is genuinely 2-parameter -- consistent with the existing Lemma C
counterexample (docs/notes/h13_line_model.md, section 4) found independently
here via a different (count-based, not inv-based) formulation of the same
obstruction.  Recorded as a negative result: do not retry a single-axis
argument for H13-I.

concordant(pi, b) = #{(j, k) : j < k, (pi[j] - b) mod n < (pi[k] - b) mod n}.
"""
import argparse
import itertools
import json
import time
from pathlib import Path


def max_b_concordant(pi, n):
    best = -1
    best_bs = []
    for b in range(n):
        cnt = 0
        for j in range(n):
            rj = (pi[j] - b) % n
            for k in range(j + 1, n):
                rk = (pi[k] - b) % n
                if rj < rk:
                    cnt += 1
        if cnt > best:
            best = cnt
            best_bs = [b]
        elif cnt == best:
            best_bs.append(b)
    return best, best_bs


def run(nmax, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for n in range(2, nmax + 1):
        t0 = time.time()
        target = (n * n) // 4
        worst_val = None
        worst_pi = None
        checked = 0
        for pi in itertools.permutations(range(n)):
            v, _ = max_b_concordant(pi, n)
            checked += 1
            if worst_val is None or v < worst_val:
                worst_val = v
                worst_pi = pi
        dt = time.time() - t0
        ok = worst_val >= target
        results.append(
            dict(n=n, target=target, min_over_pi_of_max_b=worst_val,
                 worst_pi=list(worst_pi), checked=checked, seconds=round(dt, 3),
                 holds=ok)
        )
        print(f"n={n} target={target} min-max={worst_val} "
              f"worst={worst_pi} checked={checked} {'OK' if ok else 'FAILS'} "
              f"({dt:.2f}s)")
    with open(out_dir / "report.json", "w") as f:
        json.dump(results, f, indent=2)
    return results


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--out", type=str, default="data/runs/h13i_axis_lemma")
    args = ap.parse_args()
    run(args.nmax, Path(args.out))
