"""Loss map of the N1 construction (PLAN section 4.1).

For every input pi the following chain of quantities is computed
(all exact; d(pi) from certified tables where available):

    Q_avg  = mean_c [2F_c - S_c + H(c)]          (proof side of N1 (3.4))
    Q_min  = min_c  [2F_c - S_c + H(c)]          (= min_c len of the raw word)
    Q_R    = min_c  [2F_c - S_c + R_c]           (R_c: shortest walk on Z_n
                                                  from 0 to c through all
                                                  carriers of shift c)
    Q_RJ   = min_c  reduced_len of the word assembled with the carrier-only
             route (free reduction at the junctions)
    d      = d(pi)

Losses (PLAN section 4.1, in this order):
    shift    = Q_avg - Q_min      (minimum over c versus the mean)
    route    = Q_min - Q_R        (H(c) versus the carrier-only route)
    junction = Q_R  - Q_RJ        (free cancellations at the junctions)
    residual = Q_RJ - d           (what none of the three explains)

Also recorded: Q_J = min_c reduced_len with the original route (junctions
alone), and the proof-side slacks s - rhs(7.3), s - A(x), s - B(x),
U_n - Q_avg.  The carrier-only route and the reduced words are variants
of the construction (AGENTS.md rule 16): they are measured here, no lemma
of N1 is claimed for them.

Modes:
  --perms N        all permutations 4 <= n <= N with certified tables
  --families N     reflections (all h), sigma_n rotations, C23 family, 4<=n<=N
  --sampled        n = 20, 50, 100, seed 20260912, losses versus B_n only
  --selftest       carrier-only route versus brute force, n <= 8
  --h8 N           hypothesis H8 only: Q_R - B_n and Q_RJ - B_n over all
                   permutations 4 <= n <= N (words not executed; n = 9 takes
                   minutes)
"""

import argparse
import itertools
import json
import os
import random
import sys
import time
from collections import Counter, deque
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "oracle"))
sys.path.insert(0, os.path.join(HERE, "..", "constructions"))
sys.path.insert(0, os.path.join(HERE, "..", "exact"))
sys.path.insert(0, os.path.join(HERE, "..", "bounds"))
sys.path.insert(0, os.path.join(HERE, "..", "checks"))
from moves import apply_word, freely_reduce, word_stats, delta, sigma, \
    CORE_VERSION  # noqa: E402
import strict_upper as SU  # noqa: E402
from bfs import rank_perm, factorials  # noqa: E402
from known import strict_upper_bound, target_diameter  # noqa: E402
from check_C16 import load_table, rhs_71, A_of, B_of  # noqa: E402,F401

EXPERIMENT_VERSION = "loss_map-1.0 (strict_upper-1.0, oracle-1.0)"


# ---------------------------------------------------------- carrier route ---

def carrier_route(n, c, carriers):
    """Shortest walk on the cycle Z_n from head position 0 to c that visits
    every position in `carriers`.  Returns a list of (move, position after
    move); L moves the head to pos+1, R to pos-1 (as in strict_upper.route).

    The walk covers a cyclic interval [a, b] (forward direction) containing
    0, c and all carriers; the interval is the complement of one gap between
    cyclically consecutive required points.  Within the interval the walk
    goes to one end, then to the other end, then to c."""
    c %= n
    req = sorted(set(carriers) | {0, c})
    best = None
    k = len(req)
    for i in range(k):
        a = req[(i + 1) % k]          # interval starts after the skipped gap
        b = req[i]
        L = (b - a) % n               # interval length (0 if k == 1)
        s0 = (0 - a) % n
        s1 = (c - a) % n
        cost1 = s0 + L + (L - s1)     # back to a, forward to b, back to c
        cost2 = (L - s0) + L + s1     # forward to b, back to a, forward to c
        for cost, opt in ((cost1, 1), (cost2, 2)):
            if best is None or cost < best[0]:
                best = (cost, opt, a, b, L, s0, s1)
    # full wrap: go once around the circle and on to c (needed e.g. for
    # c = 0 with every position required: n moves, not 2(n-1))
    wrap = n + delta(0, c, n)
    if wrap < best[0]:
        best = (wrap, 3, None, None, None, None, None)
    cost, opt, a, b, L, s0, s1 = best
    out, pos = [], 0

    def go(steps, direction):
        nonlocal pos
        for _ in range(steps):
            pos = (pos + direction) % n
            out.append(("L" if direction == 1 else "R", pos))

    if opt == 3:
        direction = 1 if c <= n - c else -1
        go(n + delta(0, c, n), direction)
    elif opt == 1:
        go(s0, -1)
        go(L, 1)
        go(L - s1, -1)
    else:
        go(L - s0, 1)
        go(L, -1)
        go(s1, 1)
    assert pos == c and len(out) == cost
    visited = {p for _, p in out} | {0}
    assert set(carriers) <= visited
    return out


def carrier_route_bruteforce(n, c, carriers):
    """BFS over (position, visited-set) — reference for carrier_route."""
    c %= n
    need = frozenset(carriers)
    start = (0, frozenset([0]) & need)
    dist = {start: 0}
    dq = deque([start])
    while dq:
        pos, vis = dq.popleft()
        if pos == c and vis == need:
            return dist[(pos, vis)]
        for step in (1, -1):
            q = (pos + step) % n
            nv = vis | ({q} & need)
            if (q, nv) not in dist:
                dist[(q, nv)] = dist[(pos, vis)] + 1
                dq.append((q, nv))
    raise AssertionError("unreachable")


def selftest(max_n=8):
    t0 = time.time()
    count = 0
    for n in range(2, max_n + 1):
        for c in range(n):
            for r in range(n + 1):
                for carriers in itertools.combinations(range(n), r):
                    w = carrier_route(n, c, carriers)
                    assert len(w) == carrier_route_bruteforce(n, c, carriers)
                    count += 1
    print(f"selftest: carrier_route == brute force on {count} cases, "
          f"2<=n<={max_n}, {time.time()-t0:.1f}s")


# ---------------------------------------------------------- word assembly ---

def assemble(pi, shift, walk):
    """Word of N1 section 3 for shift c with an arbitrary external walk
    (list of (move, position)); local words applied at the first visit of
    their carrier.  Returns (word, walk length)."""
    by_carrier = {lw["carrier"]: lw for lw in shift["locals"]}
    parts, used = [], set()
    if 0 in by_carrier:
        parts.append(by_carrier[0]["word"])
        used.add(0)
    for mv, pos in walk:
        parts.append(mv)
        if pos in by_carrier and pos not in used:
            parts.append(by_carrier[pos]["word"])
            used.add(pos)
    assert used == set(by_carrier), "walk missed a carrier"
    return "".join(parts), len(walk)


# --------------------------------------------------------------- per input ---

def analyse(pi, d=None, execute=True):
    n = len(pi)
    shifts = SU.all_shifts(pi)
    sm = SU.perm_summary(pi, shifts)
    ident = tuple(range(n))
    per_c = []
    for s in shifts:
        carriers = [lw["carrier"] for lw in s["locals"]]
        walk = carrier_route(n, s["c"], carriers)
        wordR, R = assemble(pi, s, walk)
        assert len(wordR) == s["local_len"] + R
        redR = freely_reduce(wordR)
        if execute:
            assert apply_word(pi, wordR) == ident, "carrier-route word sorts"
            assert apply_word(pi, redR) == ident
        lR, xR, rotR = word_stats(wordR)
        rlR, rxR, rrotR = word_stats(redR)
        per_c.append({
            "c": s["c"], "F_c": s["F_c"], "S_c": s["S_c"], "sumM": s["sumM"],
            "E_c": s["E_c"], "theta_c": s["theta_c"], "H": s["H"],
            "route_len": s["route_len"], "local_len": s["local_len"],
            "K_c": len(s["cycles"]), "carriers": carriers,
            "len": s["len"], "N_X": s["N_X"], "N_rot": s["N_rot"],
            "reduced_len": s["reduced_len"], "reduced_N_X": s["reduced_N_X"],
            "reduced_N_rot": s["reduced_N_rot"],
            "R_c": R, "bound_R": s["local_bound"] + R,
            "lenR": lR, "N_X_R": xR, "N_rot_R": rotR,
            "reducedR_len": rlR, "reducedR_N_X": rxR, "reducedR_N_rot": rrotR,
        })
    Q_avg = sm["avg_bound"]
    Q_min = sm["min_bound"]
    Q_J = sm["min_reduced_len"]
    Q_R = min(p["bound_R"] for p in per_c)
    Q_RJ = min(p["reducedR_len"] for p in per_c)
    ref = d if d is not None else target_diameter(n)
    s_ = sm["s"]
    x = sm["x"]
    rhs73 = (Fraction(4, 3) * sm["mu"] + Fraction(sm["P"], 3 * n)
             + Fraction(7, 9) * x + Fraction(4, 3) * sm["theta_bar"])
    out = {
        "pi": list(pi), "n": n, "d": d, "B_n": target_diameter(n),
        "reference": "d(pi)" if d is not None else "B_n",
        "Q_avg": str(Q_avg), "Q_min": Q_min, "Q_J": Q_J, "Q_R": Q_R,
        "Q_RJ": Q_RJ,
        "loss_shift": str(Q_avg - Q_min), "loss_route": Q_min - Q_R,
        "loss_junction": Q_R - Q_RJ, "loss_residual": Q_RJ - ref,
        "loss_junction_only": Q_min - Q_J,
        "loss_total_min": Q_min - ref,
        "s": str(s_), "x": str(x), "V": sm["V"], "T_min": sm["T_min"],
        "slack_73": str(s_ - rhs73), "slack_81": str(s_ - A_of(x, n)),
        "slack_82": str(s_ - B_of(x, n)),
        "slack_U": str(strict_upper_bound(n) - Q_avg),
        "argmin_c_bound": [p["c"] for p in per_c if p["bound_R"] == Q_R],
        "per_c": per_c,
    }
    return out


# --------------------------------------------------------------- aggregate ---

class Agg:
    KEYS = ("loss_shift", "loss_route", "loss_junction", "loss_residual",
            "loss_junction_only", "loss_total_min", "slack_U", "slack_73")

    def __init__(self, n, keep=6):
        self.n = n
        self.count = 0
        self.sums = {k: Fraction(0) for k in self.KEYS}
        self.maxs = {k: None for k in self.KEYS}
        self.hist = {k: Counter() for k in ("loss_shift_floor", "loss_route",
                                            "loss_junction", "loss_residual",
                                            "loss_total_min",
                                            "final_minus_ref")}
        self.keep = keep
        self.worst_min = []      # by Q_min - ref
        self.worst_final = []    # by Q_RJ - ref
        self.max_Q = {"Q_min": 0, "Q_R": 0, "Q_RJ": 0, "Q_J": 0}

    def add(self, r):
        self.count += 1
        for k in self.KEYS:
            v = Fraction(r[k])
            self.sums[k] += v
            if self.maxs[k] is None or v > self.maxs[k]:
                self.maxs[k] = v
        self.hist["loss_shift_floor"][int(Fraction(r["loss_shift"]))] += 1
        for k in ("loss_route", "loss_junction", "loss_residual",
                  "loss_total_min"):
            self.hist[k][r[k]] += 1
        self.hist["final_minus_ref"][r["Q_RJ"] - (r["d"] if r["d"] is not None
                                                  else r["B_n"])] += 1
        for k in self.max_Q:
            self.max_Q[k] = max(self.max_Q[k], r[k])
        slim = {k: v for k, v in r.items() if k != "per_c"}
        slim["per_c"] = [{k2: v2 for k2, v2 in p.items()} for p in r["per_c"]]
        self._push(self.worst_min, r["loss_total_min"], slim)
        self._push(self.worst_final, r["loss_residual"], slim)

    def _push(self, lst, key, item):
        lst.append((key, item))
        lst.sort(key=lambda t: (-t[0], t[1]["pi"]))
        del lst[self.keep:]

    def summary(self):
        cnt = self.count
        return {
            "n": self.n, "count": cnt,
            "mean": {k: str(self.sums[k] / cnt) for k in self.KEYS},
            "max": {k: str(self.maxs[k]) for k in self.KEYS},
            "hist": {k: dict(sorted(v.items())) for k, v in self.hist.items()},
            "max_Q": self.max_Q, "B_n": target_diameter(self.n),
            "U_n": strict_upper_bound(self.n),
            "worst_by_min_c_len": [it for _, it in self.worst_min],
            "worst_by_final": [it for _, it in self.worst_final],
        }


def fmt_summary(sm):
    n = sm["n"]
    m = sm["mean"]
    mx = sm["max"]
    ln = [f"n={n}: {sm['count']} inputs, B_n={sm['B_n']}, U_n={sm['U_n']}, "
          f"max Q_min={sm['max_Q']['Q_min']}, max Q_R={sm['max_Q']['Q_R']}, "
          f"max Q_RJ={sm['max_Q']['Q_RJ']}, max Q_J={sm['max_Q']['Q_J']}"]
    for k in Agg.KEYS:
        ln.append(f"  {k:20s} mean={float(Fraction(m[k])):7.3f}  "
                  f"max={mx[k]}")
    for k in ("loss_route", "loss_junction", "loss_residual",
              "loss_total_min", "final_minus_ref"):
        ln.append(f"  hist {k}: {sm['hist'][k]}")
    return "\n".join(ln)


# -------------------------------------------------------------------- modes ---

def run_perms(max_n, report):
    out = {}
    for n in range(4, max_n + 1):
        t0 = time.time()
        table = load_table(n)
        fact = factorials(n)
        agg = Agg(n)
        for pi in itertools.permutations(range(n)):
            d = table[rank_perm(pi, fact)] if table is not None else None
            agg.add(analyse(pi, d))
        sm = agg.summary()
        sm["seconds"] = round(time.time() - t0, 1)
        sm["coverage"] = "all permutations, all shifts, words executed"
        sm["table"] = table is not None
        out[n] = sm
        print(fmt_summary(sm) + f"\n  {sm['seconds']}s")
    report["perms"] = out


def c23_perm(b):
    n = 8 * b
    tau = [0] * n
    for k in range(8):
        for j in range(b):
            tau[b * k + j] = b * ((3 * k) % 8) + j
    return tuple((-t) % n for t in tau)


def run_families(max_n, report):
    out = {}
    for n in range(4, max_n + 1):
        table = load_table(n)
        fact = factorials(n)
        fam = {}
        # reflections pi(i) = h - i, all h (h = 1 is sigma_n)
        rows = []
        for h in range(n):
            pi = tuple((h - i) % n for i in range(n))
            d = table[rank_perm(pi, fact)] if table is not None else None
            r = analyse(pi, d, execute=n <= 40)
            rows.append({"h": h, "d": d, "Q_min": r["Q_min"], "Q_R": r["Q_R"],
                         "Q_RJ": r["Q_RJ"], "Q_J": r["Q_J"],
                         "loss_shift": r["loss_shift"],
                         "loss_route": r["loss_route"],
                         "loss_junction": r["loss_junction"],
                         "loss_residual": r["loss_residual"]})
        fam["reflections"] = rows
        # rotations of sigma_n: L^k applied to sigma_n and sigma_n L^k
        rows = []
        sg = sigma(n)
        for k in range(n):
            for kind in ("left", "right"):
                if kind == "left":
                    pi = tuple(sg[(i + k) % n] for i in range(n))
                else:
                    pi = tuple((sg[i] + k) % n for i in range(n))
                d = table[rank_perm(pi, fact)] if table is not None else None
                r = analyse(pi, d, execute=n <= 40)
                rows.append({"kind": kind, "k": k, "d": d,
                             "Q_min": r["Q_min"], "Q_R": r["Q_R"],
                             "Q_RJ": r["Q_RJ"],
                             "loss_residual": r["loss_residual"]})
        fam["sigma_rotations"] = rows
        if n % 8 == 0:
            pi = c23_perm(n // 8)
            d = table[rank_perm(pi, fact)] if table is not None else None
            r = analyse(pi, d, execute=n <= 40)
            fam["c23"] = {k: v for k, v in r.items() if k != "per_c"}
        out[n] = fam
        refl = fam["reflections"]
        print(f"families n={n} (table={table is not None}): reflections "
              f"Q_min-ref: {[r['Q_min'] - (r['d'] if r['d'] is not None else target_diameter(n)) for r in refl]}; "
              f"Q_RJ-ref: {[r['loss_residual'] for r in refl]}; "
              f"sigma rotations max residual "
              f"{max(r['loss_residual'] for r in fam['sigma_rotations'])}"
              + (f"; C23 Q_min-B_n={fam['c23']['loss_total_min']} "
                 f"Q_RJ-B_n={fam['c23']['loss_residual']}" if "c23" in fam
                 else ""))
    report["families"] = out


def run_sampled(report, seed=20260912, cases=40, orders=(20, 50, 100)):
    rng = random.Random(seed)
    out = {}
    for n in orders:
        t0 = time.time()
        agg = Agg(n, keep=3)
        for idx in range(cases):
            if idx % 2 == 0:
                pi = tuple(rng.sample(range(n), n))
            else:
                h = rng.randrange(n)
                pi = [(h - i) % n for i in range(n)]
                for _ in range(rng.randrange(1, 4)):
                    i, j = rng.randrange(n), rng.randrange(n)
                    pi[i], pi[j] = pi[j], pi[i]
                pi = tuple(pi)
            agg.add(analyse(pi, None, execute=True))
        sm = agg.summary()
        sm["seconds"] = round(time.time() - t0, 1)
        sm["coverage"] = f"SAMPLED {cases} (random + perturbed reflections), seed {seed}; reference B_n, not d(pi)"
        for key in ("worst_by_min_c_len", "worst_by_final"):
            for it in sm[key]:
                for p in it["per_c"]:
                    p.pop("carriers", None)
        out[n] = sm
        print(fmt_summary(sm) + f"\n  {sm['seconds']}s")
    report["sampled"] = out


def run_h8(max_n, report):
    """H8: min_c [2F_c - S_c + R_c] <= B_n + 1 (and <= B_n after free
    reduction) over all permutations; exceptions above B_n are listed."""
    out = {}
    for n in range(4, max_n + 1):
        t0 = time.time()
        B = target_diameter(n)
        hist_R, hist_RJ, exc = Counter(), Counter(), []
        count = 0
        for pi in itertools.permutations(range(n)):
            r = analyse(pi, None, execute=False)
            count += 1
            hist_R[r["Q_R"] - B] += 1
            hist_RJ[r["Q_RJ"] - B] += 1
            if r["Q_R"] > B:
                exc.append({"pi": list(pi), "Q_R": r["Q_R"],
                            "Q_RJ": r["Q_RJ"]})
        out[n] = {"count": count, "B_n": B,
                  "max_Q_R_minus_B": max(hist_R),
                  "max_Q_RJ_minus_B": max(hist_RJ),
                  "hist_Q_R_minus_B_top": {k: hist_R[k] for k in
                                           sorted(hist_R)[-4:]},
                  "hist_Q_RJ_minus_B_top": {k: hist_RJ[k] for k in
                                            sorted(hist_RJ)[-4:]},
                  "exceptions_Q_R_above_B": exc,
                  "seconds": round(time.time() - t0, 1)}
        print(f"H8 n={n}: {count} perms, max Q_R-B_n={max(hist_R)}, "
              f"max Q_RJ-B_n={max(hist_RJ)}, top hist Q_R-B_n "
              f"{out[n]['hist_Q_R_minus_B_top']}, exceptions {len(exc)}, "
              f"{out[n]['seconds']}s")
        for e in exc:
            print("   ", e)
    report["h8"] = out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perms", type=int)
    ap.add_argument("--families", type=int)
    ap.add_argument("--sampled", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--h8", type=int)
    ap.add_argument("--report")
    args = ap.parse_args()
    assert apply_word((2, 0, 1), "XL") == (2, 1, 0)
    report = {"core_version": CORE_VERSION,
              "construction_version": SU.CONSTRUCTION_VERSION,
              "experiment_version": EXPERIMENT_VERSION,
              "command": " ".join(sys.argv)}
    if args.selftest:
        selftest()
    if args.perms:
        run_perms(args.perms, report)
    if args.families:
        run_families(args.families, report)
    if args.sampled:
        run_sampled(report)
    if args.h8:
        run_h8(args.h8, report)
    if args.report:
        os.makedirs(os.path.dirname(args.report), exist_ok=True)
        with open(args.report, "w") as f:
            json.dump(report, f, indent=1)
    print("DONE")


if __name__ == "__main__":
    main()
