"""h3_glue_probe-1.0 -- first measurement for the H3 induction (PLAN §4.3, п. 3).

Scheme under test ("gluing"): to solve the n-problem, first bring the state to
one where the values n-1 and 0 sit in cyclically adjacent cells in that order
(as they do in the identity, at positions n-1 and 0); then treat that pair as
one super-cell and solve the (n-1)-problem.

Two costs decide whether the scheme can fit the zero-slack recurrence
`T_n <= T_{n-1} + (n - 1)`:

  S(n) = max_pi  min_{sigma in G_n} dist(pi, sigma)      -- cost of gluing,
  E(n) = max_{sigma in G_n} [ d_n(sigma) - d_{n-1}(c(sigma)) ]
                                                          -- cost of imitation,

where G_n is the set of glued states and `c(sigma)` deletes the entry n-1 from
the word (leaving a word of length n-1 over the values 0..n-2).  S(n) + E(n)
is only a crude bound (the worst gluing target and the worst imitation state
need not be the same); the decisive quantity is the joint optimum

  U(n) = max_pi min_{sigma in G_n} [ dist(pi, sigma) + d_{n-1}(c(sigma)) ],

which turns out to be vacuous (U(pi) <= d_n(pi) for every pi, because the
imitation is assumed free), and the quantity that really binds the recurrence

  R(n) = max_pi min_{sigma in G_n} [ dist(pi, sigma) + E(sigma) ],
         E(sigma) = d_n(sigma) - d_{n-1}(c(sigma)).

Indeed d_n(pi) <= dist(pi, sigma) + d_n(sigma) = dist(pi, sigma) + E(sigma)
+ d_{n-1}(c(sigma)) <= min_sigma[...] + T_{n-1}, so `R(n) <= n - 1` closes
`T_n <= T_{n-1} + (n - 1)` and `R(n) > n - 1` refutes this gluing as it stands.

Distances are computed by BFS with the reference moves (oracle-1.0), not read
from the tables, so the probe is self-contained; the values of D_n are printed
and can be compared with B_n by eye.

Usage: python3 experiments/h3_glue_probe.py [--nmax 8]
"""

import argparse
import itertools
import json
import os
import sys
import time
from collections import deque

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))

from moves import apply_move, identity, MOVES, CORE_VERSION  # noqa: E402

VERSION = "h3_glue_probe-1.0"


def bfs(n, sources):
    """Multi-source BFS over S_n with the LRX generators."""
    dist = {}
    dq = deque()
    for s in sources:
        dist[s] = 0
        dq.append(s)
    while dq:
        p = dq.popleft()
        d = dist[p]
        for m in MOVES:
            q = apply_move(p, m)
            if q not in dist:
                dist[q] = d + 1
                dq.append(q)
    return dist


def glued(n):
    """States where n-1 is immediately followed (cyclically) by 0."""
    out = []
    for p in itertools.permutations(range(n)):
        i = p.index(n - 1)
        if p[(i + 1) % n] == 0:
            out.append(p)
    return out


def contract(state, n):
    """Delete the entry n-1; the result is a word of length n-1 over 0..n-2."""
    return tuple(v for v in state if v != n - 1)


def glued_any(n):
    """States where some value v is immediately followed (cyclically) by
    v+1 mod n -- the weaker gluing: any consecutive pair, not only (n-1, 0)."""
    out = []
    for p in itertools.permutations(range(n)):
        for i in range(n):
            if p[(i + 1) % n] == (p[i] + 1) % n:
                out.append(p)
                break
    return out


def contract_any(state, n):
    """For a state of glued_any: delete the value v+1 of the first glued pair
    and relabel the survivors by rank.  The identity maps to the identity."""
    for i in range(n):
        if state[(i + 1) % n] == (state[i] + 1) % n:
            drop = (state[i] + 1) % n
            return tuple(v - 1 if v > drop else v for v in state if v != drop)
    raise ValueError("not glued")


def joint_optimum(n, G, weight):
    """For every state pi: min over sigma in G of dist(pi, sigma) + weight(sigma).
    Unit edge weights, integer initial potentials -> bucket-queue Dijkstra."""
    best = {}
    buckets = {}
    for s in G:
        w = weight(s)
        if best.get(s, 10 ** 9) > w:
            best[s] = w
            buckets.setdefault(w, []).append(s)
    k = min(buckets)
    kmax = max(buckets)
    while k <= kmax:
        for p in buckets.get(k, ()):
            if best[p] != k:
                continue
            for m in MOVES:
                q = apply_move(p, m)
                if best.get(q, 10 ** 9) > k + 1:
                    best[q] = k + 1
                    buckets.setdefault(k + 1, []).append(q)
                    kmax = max(kmax, k + 1)
        buckets.pop(k, None)
        k += 1
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    out = args.out or os.path.join(ROOT, "data", "runs", "h3_glue_probe")
    os.makedirs(out, exist_ok=True)

    rows = []
    prev = None          # distance table of n-1
    for n in range(3, args.nmax + 1):
        t0 = time.time()
        d = bfs(n, [identity(n)])
        D = max(d.values())
        row = {"n": n, "D_n": D, "B_n": n * (n - 1) // 2}
        if prev is not None:
            G = glued(n)
            g = bfs(n, G)
            S = max(g.values())
            argS = max(g, key=lambda p: (g[p], tuple(-x for x in p)))
            E = -10 ** 9
            argE = None
            for s in G:
                e = d[s] - prev[contract(s, n)]
                if e > E:
                    E, argE = e, s
            U = joint_optimum(n, G, lambda s: prev[contract(s, n)])
            Umax = max(U.values())
            argU = max(U, key=lambda p: (U[p], tuple(-x for x in p)))
            loss = max(U[p] - d[p] for p in U)
            argL = max(U, key=lambda p: (U[p] - d[p], tuple(-x for x in p)))
            Rm = joint_optimum(n, G, lambda s: d[s] - prev[contract(s, n)])
            Rmax = max(Rm.values())
            argR = max(Rm, key=lambda p: (Rm[p], tuple(-x for x in p)))
            G2 = glued_any(n)
            g2 = bfs(n, G2)
            S2 = max(g2.values())
            R2m = joint_optimum(n, G2, lambda s: d[s] - prev[contract_any(s, n)])
            R2 = max(R2m.values())
            argR2 = max(R2m, key=lambda p: (R2m[p], tuple(-x for x in p)))
            # refined form: does every state have a shortest word through G2?
            P2m = joint_optimum(n, G2, lambda s: d[s])
            P2 = max(P2m[p] - d[p] for p in P2m)
            argP2 = max(P2m, key=lambda p: (P2m[p] - d[p], tuple(-x for x in p)))
            Q2 = max(P2m.values())
            row.update({"detour_any": P2, "arg_detour_any": list(argP2),
                        "max_through_any": Q2,
                        "through_any_fits": Q2 <= n * (n - 1) // 2})
            row.update({"any_glued_states": len(G2), "S_any": S2, "R_any": R2,
                        "arg_R_any": list(argR2), "R_any_fits": R2 <= n - 1,
                        "R_any_minus_budget": R2 - (n - 1)})
            row.update({"glued_states": len(G), "S": S, "arg_S": list(argS),
                        "E": E, "arg_E": list(argE), "S_plus_E": S + E,
                        "budget": n - 1, "crude_fits": S + E <= n - 1,
                        "U": Umax, "arg_U": list(argU),
                        "U_minus_Bn": Umax - n * (n - 1) // 2,
                        "U_fits": Umax <= n * (n - 1) // 2,
                        "max_loss_vs_d": loss, "arg_loss": list(argL),
                        "R": Rmax, "arg_R": list(argR),
                        "R_fits": Rmax <= n - 1, "R_minus_budget": Rmax - (n - 1)})
        row["seconds"] = round(time.time() - t0, 2)
        rows.append(row)
        print(row)
        prev = d

    res = {"version": VERSION, "core": CORE_VERSION, "rows": rows}
    with open(os.path.join(out, "probe.json"), "w") as f:
        json.dump(res, f, indent=1)
    print("written:", os.path.join(out, "probe.json"))


if __name__ == "__main__":
    main()
