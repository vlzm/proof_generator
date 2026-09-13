"""Ad hoc SAMPLED check (not a permanent module): does Q_R stay <= B_n on
sigma_n, rev_n, the C23 family and a few random permutations at n well
beyond the exhaustive families range (4<=n<=12) of experiments/loss_map.py?
Uses only existing, audited code (experiments/loss_map.py: analyze,
constructions/strict_upper.py, oracle/moves.py) -- no new logic.

Run from the repository root:
    PYTHONPATH= python3 data/runs/h10_reflection_route/other_families_check.py
Seeds: random.Random(42+n) for n<=40, random.Random(1000+n) for n>40.
"""
import os
import random
import sys
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(ROOT, "..", "..", "..")
for sub in ("oracle", "exact", "bounds", "constructions", "experiments"):
    sys.path.insert(0, os.path.join(ROOT, sub))
from moves import sigma, rev  # noqa: E402
import known  # noqa: E402
import loss_map as lm  # noqa: E402


def check(name, pi, n):
    t0 = time.time()
    r = lm.analyze(tuple(pi))
    bn = known.target_diameter(n)
    print(f"n={n:4d} {name:16s} Q_R-B={r['Q_R']-bn:+d} Q_min-B={r['Q_min']-bn:+d} time={time.time()-t0:.2f}s")


def c23_family(n):
    b = n // 8
    tau = [0] * n
    for k in range(8):
        for j in range(b):
            tau[b * k + j] = b * ((3 * k) % 8) + j
    return tuple((-tau[i]) % n for i in range(n))


def run(ns, seed_base, n_random):
    for n in ns:
        check("sigma_n", sigma(n), n)
        check("rev_n", rev(n), n)
        if n % 8 == 0:
            check("C23_family", c23_family(n), n)
        rng = random.Random(seed_base + n)
        for t in range(n_random):
            p = list(range(n))
            rng.shuffle(p)
            check(f"random{t}", p, n)


if __name__ == "__main__":
    run((16, 24, 32, 40), 42, 3)
    run((64, 128, 200), 1000, 2)
