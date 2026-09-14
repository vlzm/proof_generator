"""H13-I, session 9: three cheap heuristics for choosing the double cut (a, b)
that all FAIL to reach the exhaustive bound floor((n-1)^2/4) for some n <= 8,
even though the true minimum I(pi) = min_{a,b} f(a,b) never exceeds it
(experiments/torus_cut.c, checks/check_C37.py). Recorded so future sessions
do not re-try the same shortcuts (see docs/notes/h13_line_model.md §7).

  H1 "fix a = 0, vary b only": min_b f(0, b). Fails already at n = 7.
  H2 "greedy two-step": a* = argmin_a f(a, 0) (best single cut, b = 0), then
     min_b f(a*, b). Fails at n = 7 too (a strictly smaller counterexample
     than H1's n = 7 one is not implied and not searched for here).
  H3 "canonical b = pi(a)" (origin at a literal point of the permutation
     matrix): min_a f(a, pi(a)). Fails badly (grows like C(n,2), i.e. no
     better than not cutting at all) already at n = 4.

Usage: python3 experiments/h13i_heuristics.py [--nmax 8]
Output: data/runs/h13i_heuristics/report.json, report.md.
"""

import argparse
import itertools
import json
import os
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
VERSION = "h13i_heuristics-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h13i_heuristics")


def brute_inv(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def cut_inv(pi, a, b):
    n = len(pi)
    return brute_inv([(pi[(a + j) % n] - b) % n for j in range(n)])


def h1_fixed_a0(pi):
    n = len(pi)
    return min(cut_inv(pi, 0, b) for b in range(n))


def h2_greedy_two_step(pi):
    n = len(pi)
    astar = min(range(n), key=lambda a: cut_inv(pi, a, 0))
    return min(cut_inv(pi, astar, b) for b in range(n))


def h3_canonical_b(pi):
    n = len(pi)
    return min(cut_inv(pi, a, pi[a]) for a in range(n))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    lines = []

    def log(msg):
        print(msg, flush=True)
        lines.append(msg)

    log(f"== {VERSION} args={vars(args)}")
    heuristics = {"H1_fixed_a0": h1_fixed_a0, "H2_greedy_two_step": h2_greedy_two_step,
                  "H3_canonical_b": h3_canonical_b}
    rows = []
    for n in range(4, args.nmax + 1):
        t0 = time.time()
        bound = (n - 1) ** 2 // 4
        worst = {name: (-1, None) for name in heuristics}
        for perm in itertools.permutations(range(n)):
            for name, fn in heuristics.items():
                v = fn(list(perm))
                if v > worst[name][0]:
                    worst[name] = (v, perm)
        row = {"n": n, "bound": bound, "seconds": round(time.time() - t0, 2)}
        for name, (v, perm) in worst.items():
            row[name] = {"worst": v, "witness": list(perm), "exceeds_bound": v > bound}
            log(f"n={n} {name}: worst over all pi = {v} (bound {bound}) "
                f"{'EXCEEDS' if v > bound else 'within bound'}; witness {perm}")
        rows.append(row)
        log(f"n={n}: {time.time() - t0:.2f} s")
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump({"version": VERSION, "args": vars(args), "rows": rows}, f, indent=1)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# h13i_heuristics — три неудачные эвристики выбора (a, b) для H13-I\n\n")
        f.write(f"Команда: `python3 experiments/h13i_heuristics.py --nmax {args.nmax}`. "
                f"Версия: {VERSION}.\n\n```text\n")
        f.write("\n".join(lines) + "\n```\n")


if __name__ == "__main__":
    main()
