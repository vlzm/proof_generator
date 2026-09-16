"""H3 induction probe (session 9): does the natural "freeze the top element"
embedding S_{n-1} -> S_n realize the recurrence T_n <= T_{n-1} + (n-1)?

Embedding (PLAN.md §4.3 item 3): for pi in S_n, let k = pi.index(n-1) (the
position of the largest value).  Rotate pi by the shorter direction so that
value n-1 lands at position n-1: rot_cost = delta_n(k, n-1) <= floor(n/2).
The rotated array p2 then has p2[n-1] = n-1 and p2[0:n-1] = q, a genuine
element of S_{n-1} (the remaining n-1 values read starting right after where
n-1 used to sit on the circle -- this part is exact, not an approximation).

Candidate construction ("naive verbatim lift"): take Wsort, a shortest word
sorting q to id_{n-1} in the (n-1)-generators, and apply the SAME letters
(L, R, X) directly to p2 as n-generators.  Since L/R on the big array always
move ALL n entries (PLAN.md's flagged difficulty: "a fixed element cannot be
silently excluded"), value n-1 drifts around the circle one slot per L/R and
can collide with the small word's X operations (which always act on the
current head, positions 0/1).  This script measures how far short that
naive construction falls: rot_cost + len(Wsort) + (distance from the actual
resulting array back to id_n, looked up by exact BFS).

Finding (see docs/notes/h3_induction.md): the naive lift does NOT stay within
budget in general.  On sigma_n it happens to equal B_n at n = 5, 6 but
overshoots to B_n + 3 (n = 7) and B_n + 2 (n = 8); i.e. this specific
construction is REFUTED as a proof strategy for H3, independently of whether
H3 itself is true.  Register: CLAIMS.md C38.

Usage: python3 experiments/h3_induction_probe.py --nmax 9 [--samples 20]
Exhaustive BFS (both the (n-1)-subproblem and the full n-problem, run fresh
each time -- no dependency on data/tables/) restricts n to <= 9 in reasonable
time; --samples adds random permutations per n in addition to sigma_n.
Version h3_induction_probe-1.0.
"""

import argparse
import itertools
import json
import os
import random
import sys
import time
from collections import deque

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
from moves import apply_move, apply_word, delta, identity, sigma, MOVES, INVERSE, CORE_VERSION  # noqa: E402

VERSION = "h3_induction_probe-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h3_induction_probe")


def full_bfs(n):
    """Exact distances from id_n over the whole S_n (fresh BFS, small n only)."""
    start = identity(n)
    dist = {start: 0}
    dq = deque([start])
    while dq:
        s = dq.popleft()
        d = dist[s]
        for mv in MOVES:
            s2 = apply_move(s, mv)
            if s2 not in dist:
                dist[s2] = d + 1
                dq.append(s2)
    return dist


def shortest_word_to(target, dist_and_prev_cache):
    """Shortest word id_m -> target, reusing a cached (dist, prev) BFS tree."""
    dist, prev = dist_and_prev_cache
    word = []
    s = target
    while prev[s][0] is not None:
        p, mv = prev[s]
        word.append(mv)
        s = p
    return "".join(reversed(word))


def bfs_tree(m):
    start = identity(m)
    prev = {start: (None, None)}
    dq = deque([start])
    while dq:
        s = dq.popleft()
        for mv in MOVES:
            s2 = apply_move(s, mv)
            if s2 not in prev:
                prev[s2] = (s, mv)
                dq.append(s2)
    dist = {s: 0 for s in prev}  # unused, placeholder for signature symmetry
    return dist, prev


def naive_lift(pi, small_tree):
    """Apply the construction described in the module docstring; return a dict
    of diagnostics.  small_tree is bfs_tree(n - 1) (built once per n)."""
    n = len(pi)
    k = pi.index(n - 1)
    t_left = (k - (n - 1)) % n
    t_right = (n - 1 - k) % n
    if t_left <= t_right:
        rot_cost, rot_word = t_left, "L" * t_left
    else:
        rot_cost, rot_word = t_right, "R" * t_right
    p2 = apply_word(pi, rot_word)
    assert p2[n - 1] == n - 1
    q = p2[: n - 1]
    Wq = shortest_word_to(q, small_tree)  # id_{n-1} -> q
    d_small = len(Wq)
    Wsort = "".join(INVERSE[m] for m in reversed(Wq))  # q -> id_{n-1}
    assert apply_word(q, Wsort) == identity(n - 1)
    final = apply_word(p2, Wsort)
    return {
        "rot_cost": rot_cost,
        "d_small": d_small,
        "final": final,
        "naive_total": rot_cost + d_small,
    }


def run(nmax, samples, seed):
    random.seed(seed)
    os.makedirs(OUT, exist_ok=True)
    rows = []
    for n in range(5, nmax + 1):
        t0 = time.time()
        big_dist = full_bfs(n)
        small_tree = bfs_tree(n - 1)
        Bn = n * (n - 1) // 2
        Dn = big_dist[sigma(n)]
        inputs = [("sigma_n", sigma(n))]
        pool = list(itertools.permutations(range(n)))
        for i in range(samples):
            inputs.append((f"random_{i}", random.choice(pool)))
        for label, pi in inputs:
            r = naive_lift(pi, small_tree)
            extra = big_dist[r["final"]]
            total = r["rot_cost"] + r["d_small"] + extra
            rows.append({
                "n": n, "label": label, "pi": list(pi), "d_pi": big_dist[pi],
                "rot_cost": r["rot_cost"], "d_small": r["d_small"],
                "extra_after_lift": extra, "naive_lift_total": total,
                "B_n": Bn, "over_Bn": total - Bn,
            })
        print(f"n={n} B_n={Bn} D_n={Dn} elapsed={time.time()-t0:.1f}s "
              f"(sigma_n naive_total={[r for r in rows if r['n']==n and r['label']=='sigma_n'][0]['naive_lift_total']})")
    report = {"version": VERSION, "core_version": CORE_VERSION, "rows": rows}
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(report, f, indent=1)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# h3_induction_probe report\n\n")
        f.write("n | label | d(pi) | rot_cost | d_small | extra_after_lift | naive_total | B_n | over_Bn\n")
        f.write("---|---|---|---|---|---|---|---|---\n")
        for r in rows:
            f.write(f"{r['n']} | {r['label']} | {r['d_pi']} | {r['rot_cost']} | {r['d_small']} | "
                    f"{r['extra_after_lift']} | {r['naive_lift_total']} | {r['B_n']} | {r['over_Bn']}\n")
    print("report written to", OUT)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--samples", type=int, default=10)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    run(args.nmax, args.samples, args.seed)
