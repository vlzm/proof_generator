"""Loss map of construction N1 (PLAN 4.1): where does the proved bound lose
against the actual distance d(pi)?

Version loss_map-1.1. Uses constructions/strict_upper.py (strict_upper-1.0),
core oracle-1.0, certified tables data/tables/dist_n{4..10}.bin (C2v).
All arithmetic exact (int / Fraction) except the irrational r_n slacks,
which are floats and marked NUMERICAL.

For every input pi and every shift c in Z_n the N1 word is built and executed.
Quantities (per c):
  local_len  = sum |W_C| = 2F_c - S_c              (N1 (2.6), exact by construction)
  local_d    = sum d(state_C)                       (n <= 10: exact distance of the
               state that W_C sorts; any word sorting it keeps the carrier at 0,
               so this is a true lower bound for the local job)
  H(c)       = N1 external cost;  R_c = shortest walk 0 -> carriers -> c (variant)
  len        = 2F_c - S_c + H(c) = |word|;  red = |freely reduced word|
  len_R      = 2F_c - S_c + R_c; red_R = reduced length of the carrier-route word
Chain (minima over all c):
  Q_avg = mean_c len = right side of (3.4) with exact s
  Q_min = min_c len = min_c bound (identical in N1)
  Q_R   = min_c len_R ;  Q_RJ = min_c red_R ;  Q_J = min_c red
Losses of PLAN 4.1:
  1 local     = local_len - local_d at the chosen shift (and per c)
  2 route     = H(c*) - R_{c*} at the chosen shift; chain term Q_min - Q_R
  3 shift     = Q_avg - Q_min  (mean vs min; min of bound = min of length here)
  4 s-chain   = slacks s -> (7.3) -> (8.1)/(8.2) -> (8.3) -> rounding
  5 junctions = len - red (N1 route) and len_R - red_R; N_X, N_rot both ways
  residual    = Q_RJ - d(pi) (composition of cycles, choice of cycles/words)
shift + route + junction + residual = Q_avg - d(pi).

Usage:
  python3 experiments/loss_map.py --perms 8            # all permutations 4<=n<=8
  python3 experiments/loss_map.py --families 12        # reflections, rotations, perturbations, C23
  python3 experiments/loss_map.py --samples            # n = 20, 50, 100 (no d)
Output: data/runs/strict_loss_audit/<part>.json and <part>.md
"""

import argparse
import itertools
import json
import math
import os
import random
import sys
import time
from fractions import Fraction

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
for sub in ("oracle", "exact", "bounds", "constructions"):
    sys.path.insert(0, os.path.join(ROOT, sub))

from moves import apply_word, freely_reduce, sigma, rev, delta, CORE_VERSION  # noqa: E402
from bfs import rank_perm, factorials  # noqa: E402
import known  # noqa: E402
import strict_upper as su  # noqa: E402

EXPERIMENT_VERSION = "loss_map-1.1"
OUT_DIR = os.path.join(ROOT, "data", "runs", "strict_loss_audit")

_TABLES = {}


def table(n):
    if n not in _TABLES:
        path = os.path.join(ROOT, "data", "tables", f"dist_n{n}.bin")
        if os.path.exists(path):
            with open(path, "rb") as f:
                _TABLES[n] = (f.read(), factorials(n))
        else:
            _TABLES[n] = None
    return _TABLES[n]


def dist(p):
    t = table(len(p))
    if t is None:
        return None
    return t[0][rank_perm(tuple(p), t[1])]


# ------------------------------------------------------------------ N1 sections 4-6 quantities

def half_window(a, n):
    return {(a + i) % n for i in range(1, n // 2 + 1)}


def proof_side(pi, per_c):
    """Exact s, mu, x, theta-bar, anchor T_h (N1 5.2), T_min and the slacks of the
    chain (3.4) -> (7.3) -> (8.1)/(8.2) -> (8.3) -> U_n."""
    n = len(pi)
    P = known.P(n)
    s = Fraction(sum(pc["S_c"] for pc in per_c), n)
    mu = Fraction(sum(pc["sumM"] for pc in per_c), n)
    Ebar = Fraction(sum(pc["E_c"] for pc in per_c), n)
    thetabar = Fraction(sum(pc["theta_c"] for pc in per_c), n)
    x = Ebar - 2 * thetabar                              # (4.5)
    tau = [(-pi[i]) % n for i in range(n)]
    windows = [half_window(i, n) for i in range(n)]
    e = [len({tau[j] for j in windows[i]} ^ windows[tau[i]]) for i in range(n)]
    V = sum(e)
    assert Fraction(V, n) == x
    m = min(e)

    def T(h):
        return sum(delta((pi[i] + i) % n, h, n) for i in range(n))
    T_anchor = min(T((a + pi[a]) % n) for a in range(n) if e[a] == m)   # best admissible anchor
    T_min = min(T(h) for h in range(n))
    rhs73 = Fraction(4, 3) * mu + Fraction(P, 3 * n) + Fraction(7, 9) * x + Fraction(4, 3) * thetabar
    A = Fraction(5 * P, 3 * n) + Fraction(7, 9) * x
    Bx = Fraction(2 * (n - 1), 3) + Fraction(P, 3 * n) - (Fraction(5, 9) - Fraction(4, 3 * n)) * x - Fraction(4, 3 * n) * x * x
    z = 3 * n * n - 4 * n + 1 - 4 * P
    r_n = (math.sqrt(z) - (n - 1)) / 2
    A_r = 5 * P / (3 * n) + 7 * r_n / 9
    R_n = 2 * P + n - 2 - 2 * P / (3 * n) + 2 / n - 7 * r_n / 9
    U = known.strict_upper_bound(n)
    expr34 = 2 * P + n - 2 + Fraction(P + 2, n) - s
    out = {
        "s": s, "mu": mu, "x": x, "bar_theta": thetabar, "V": V, "m": m,
        "T_anchor": T_anchor, "T_min": T_min,
        "expr_3_4": expr34, "U_n": U,
        # (3.4) -> (7.3)
        "g1_s_minus_73": s - rhs73,
        # (7.3) -> (8.1): drop 4/3 theta, mu >= P/n
        "g2a_73_minus_81": rhs73 - A,
        "g2a_mu_slack_64a": mu - Fraction(P, n),
        # (7.3) -> (8.2): drop 4/3 theta, mu >= (n-1)/2 - T_h/n, T_h <= (n-1)x + x^2
        "g2b_73_minus_82": rhs73 - Bx,
        "g2b_mu_slack_64b": mu - (Fraction(n - 1, 2) - Fraction(T_anchor, n)),
        "g2b_Th_slack_51": (n - 1) * x + x * x - T_anchor,
        "g2_73_minus_max81_82": rhs73 - max(A, Bx),
        "which_of_81_82_is_larger": "8.1" if A >= Bx else "8.2",
        # (8.1)/(8.2) -> (8.3)   NUMERICAL
        "g3_max_minus_A_r_n": float(max(A, Bx)) - A_r,
        # rounding (T)
        "g4_R_n_minus_U_n": R_n - U,
        "slack_U_minus_expr34": U - expr34,
        "x_minus_r_n": float(x) - r_n,
    }
    return out


# ------------------------------------------------------------------ per input

def analyze(pi, want_words=False):
    n = len(pi)
    d_pi = dist(pi)
    per_c = []
    for c in range(n):
        w1 = su.shift_word(pi, c, "n1")
        wR = su.shift_word(pi, c, "carrier")
        st = dict(w1["stats"])
        loc_d = 0
        loc_red = 0
        have_d = True
        for lw in w1["locals"]:
            dd = dist(lw["state"])
            if dd is None:
                have_d = False
            else:
                loc_d += dd
            loc_red += len(freely_reduce(lw["word"]))
        red = freely_reduce(w1["word"])
        redR = freely_reduce(wR["word"])
        st.update({
            "K_c": st.pop("n_cycles"),
            "local_d": loc_d if have_d else None,
            "local_red": loc_red,
            "loss_local": (st["local_len"] - loc_d) if have_d else None,
            "R_c": wR["stats"]["route_len"],
            "len_R": wR["stats"]["len"],
            "red": len(red), "red_N_X": red.count("X"),
            "red_R": len(redR), "red_R_N_X": redR.count("X"),
            "N_X_R": wR["stats"]["N_X"],
        })
        st["red_N_rot"] = st["red"] - st["red_N_X"]
        st["red_R_N_rot"] = st["red_R"] - st["red_R_N_X"]
        st["N_rot_R"] = st["len_R"] - st["N_X_R"]
        if want_words:
            st["word"] = w1["word"]
            st["word_R"] = wR["word"]
        per_c.append(st)
    lens = [pc["len"] for pc in per_c]
    Q_avg = Fraction(sum(lens), n)
    Q_min = min(lens)
    cstar = lens.index(Q_min)                      # smallest c among ties
    Q_R = min(pc["len_R"] for pc in per_c)
    Q_RJ = min(pc["red_R"] for pc in per_c)
    Q_J = min(pc["red"] for pc in per_c)
    ps = proof_side(pi, per_c)
    assert ps["expr_3_4"] == Q_avg, (pi, ps["expr_3_4"], Q_avg)
    res = {
        "pi": list(pi), "n": n, "d": d_pi,
        "Q_avg": Q_avg, "Q_min": Q_min, "argmin_c": cstar,
        "Q_R": Q_R, "Q_RJ": Q_RJ, "Q_J": Q_J,
        "N_X_at_argmin": per_c[cstar]["N_X"], "N_rot_at_argmin": per_c[cstar]["N_rot"],
        "red_at_argmin": per_c[cstar]["red"],
        "loss_shift": Q_avg - Q_min,
        "loss_route_chain": Q_min - Q_R,
        "loss_route_at_argmin": per_c[cstar]["H"] - per_c[cstar]["R_c"],
        "loss_junction_chain": Q_R - Q_RJ,
        "loss_junction_n1_at_argmin": per_c[cstar]["len"] - per_c[cstar]["red"],
        "loss_local_at_argmin": per_c[cstar]["loss_local"],
        "loss_local_mean_c": (Fraction(sum(pc["loss_local"] for pc in per_c), n)
                              if per_c[0]["loss_local"] is not None else None),
        "loss_local_min_c": (min(pc["loss_local"] for pc in per_c)
                             if per_c[0]["loss_local"] is not None else None),
        "min_c_local_d_plus_R": (min(pc["local_d"] + pc["R_c"] for pc in per_c)
                                 if per_c[0]["local_d"] is not None else None),
        # decomposition of the shift loss Q_avg - Q_min at c*: mean_c - value at c* of
        # 2F_c, -S_c, H(c) separately (mean F_c = P, mean H = n-2+(P+2)/n)
        "shift_decomp": {
            "2F": 2 * (Fraction(sum(pc["F_c"] for pc in per_c), n) - per_c[cstar]["F_c"]),
            "-S": -(ps["s"] - per_c[cstar]["S_c"]),
            "H": Fraction(sum(pc["H"] for pc in per_c), n) - per_c[cstar]["H"],
        },
        "proof": ps,
        "per_c": per_c,
    }
    if d_pi is not None:
        res["loss_residual"] = Q_RJ - d_pi
        res["loss_total"] = Q_avg - d_pi
        res["Q_min_minus_d"] = Q_min - d_pi
        assert res["loss_shift"] + res["loss_route_chain"] + res["loss_junction_chain"] + res["loss_residual"] == res["loss_total"]
    return res


# ------------------------------------------------------------------ aggregation

def hist_add(h, key):
    h[key] = h.get(key, 0) + 1


class Agg:
    FIELDS = ("loss_shift", "loss_route_chain", "loss_route_at_argmin", "loss_junction_chain",
              "loss_junction_n1_at_argmin", "loss_local_at_argmin", "loss_local_min_c",
              "loss_residual", "Q_min_minus_d", "loss_total")
    PROOF = ("g1_s_minus_73", "g2a_73_minus_81", "g2b_73_minus_82", "g2_73_minus_max81_82",
             "g2a_mu_slack_64a", "g2b_mu_slack_64b", "g2b_Th_slack_51",
             "g3_max_minus_A_r_n", "g4_R_n_minus_U_n", "slack_U_minus_expr34")

    def __init__(self, n):
        self.n = n
        self.count = 0
        self.hist = {f: {} for f in self.FIELDS}
        self.sum = {f: 0 for f in self.FIELDS + self.PROOF}
        self.max = {f: None for f in self.FIELDS + self.PROOF}
        self.min = {f: None for f in self.FIELDS + self.PROOF}
        self.zero = {f: 0 for f in self.PROOF}
        self.larger = {"8.1": 0, "8.2": 0}
        self.worst = {"Q_min_minus_d": [], "loss_residual": [], "loss_local_at_argmin": [], "Q_R": [], "Q_RJ": []}
        self.hist_QR_minus_B = {}
        self.hist_QRJ_minus_B = {}
        self.hist_Qmin_minus_B = {}
        self.mean_c_len_R_minus_B_max = None
        self.shift_decomp_sum = {"2F": 0, "-S": 0, "H": 0}
        self.shift_decomp_max = {"2F": None, "-S": None, "H": None}

    def add(self, r):
        self.count += 1
        Bn = known.target_diameter(self.n)
        for f in self.FIELDS:
            v = r.get(f)
            if v is None:
                continue
            hist_add(self.hist[f], str(v))
            self.sum[f] += v
            self.max[f] = v if self.max[f] is None else max(self.max[f], v)
            self.min[f] = v if self.min[f] is None else min(self.min[f], v)
        for f in self.PROOF:
            v = r["proof"][f]
            self.sum[f] += v
            self.max[f] = v if self.max[f] is None else max(self.max[f], v)
            self.min[f] = v if self.min[f] is None else min(self.min[f], v)
            if v == 0:
                self.zero[f] += 1
        self.larger[r["proof"]["which_of_81_82_is_larger"]] += 1
        hist_add(self.hist_QR_minus_B, str(r["Q_R"] - Bn))
        hist_add(self.hist_QRJ_minus_B, str(r["Q_RJ"] - Bn))
        hist_add(self.hist_Qmin_minus_B, str(r["Q_min"] - Bn))
        mR = Fraction(sum(pc["len_R"] for pc in r["per_c"]), self.n) - Bn
        if self.mean_c_len_R_minus_B_max is None or mR > self.mean_c_len_R_minus_B_max:
            self.mean_c_len_R_minus_B_max = mR
        for k, v in r["shift_decomp"].items():
            self.shift_decomp_sum[k] += v
            self.shift_decomp_max[k] = v if self.shift_decomp_max[k] is None else max(self.shift_decomp_max[k], v)
        for key in self.worst:
            v = r.get(key)
            if v is None:
                continue
            lst = self.worst[key]
            lst.append((v, r["pi"]))
            lst.sort(reverse=True)
            del lst[6:]

    def summary(self):
        def fr(v):
            return str(v) if isinstance(v, (int, Fraction)) else (None if v is None else round(v, 4))
        out = {"n": self.n, "count": self.count, "B_n": known.target_diameter(self.n),
               "U_n": known.strict_upper_bound(self.n)}
        for f in self.FIELDS:
            if self.max[f] is None:
                continue
            out[f] = {"mean": fr(Fraction(self.sum[f], self.count)), "max": fr(self.max[f]),
                      "min": fr(self.min[f]), "hist": dict(sorted(self.hist[f].items(), key=lambda kv: Fraction(kv[0])))}
        for f in self.PROOF:
            mean = self.sum[f] / self.count
            out[f] = {"mean": fr(mean), "max": fr(self.max[f]), "min": fr(self.min[f]), "zero_count": self.zero[f]}
        out["which_of_81_82_binds"] = self.larger
        out["hist_Q_min_minus_B_n"] = dict(sorted(self.hist_Qmin_minus_B.items(), key=lambda kv: int(kv[0])))
        out["hist_Q_R_minus_B_n"] = dict(sorted(self.hist_QR_minus_B.items(), key=lambda kv: int(kv[0])))
        out["hist_Q_RJ_minus_B_n"] = dict(sorted(self.hist_QRJ_minus_B.items(), key=lambda kv: int(kv[0])))
        out["max_over_pi_of_mean_c_len_R_minus_B_n"] = fr(self.mean_c_len_R_minus_B_max)
        out["shift_loss_decomposition_at_argmin"] = {k: {"mean": fr(Fraction(self.shift_decomp_sum[k], self.count)), "max": fr(self.shift_decomp_max[k])}
                                                     for k in self.shift_decomp_sum}
        out["worst"] = {k: [(fr(v), p) for v, p in lst] for k, lst in self.worst.items()}
        return out


def to_json(obj):
    if isinstance(obj, Fraction):
        return str(obj)
    raise TypeError(type(obj))


def rule15_table(r):
    """Markdown table of all per-shift statistics (AGENTS.md rule 15) for one input."""
    lines = [f"pi = {tuple(r['pi'])}, n = {r['n']}, d = {r['d']}, Q_avg = {r['Q_avg']}, Q_min = {r['Q_min']} (c* = {r['argmin_c']}), "
             f"Q_J = {r['Q_J']}, Q_R = {r['Q_R']}, Q_RJ = {r['Q_RJ']}; s = {r['proof']['s']}, x = {r['proof']['x']}, "
             f"T_anchor = {r['proof']['T_anchor']}, T_min = {r['proof']['T_min']}",
             "",
             "| c | K | F_c | S_c | sum M | E_c | theta | H | R | carriers | local_len | local_d | local_red | len | red | len_R | red_R | N_X | N_rot | red N_X | red N_rot |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for pc in r["per_c"]:
        lines.append(f"| {pc['c']} | {pc['K_c']} | {pc['F_c']} | {pc['S_c']} | {pc['sumM']} | {pc['E_c']} | {pc['theta_c']} | {pc['H']} | {pc['R_c']} | "
                     f"{','.join(map(str, pc['carriers']))} | {pc['local_len']} | {pc['local_d']} | {pc['local_red']} | {pc['len']} | {pc['red']} | "
                     f"{pc['len_R']} | {pc['red_R']} | {pc['N_X']} | {pc['N_rot']} | {pc['red_N_X']} | {pc['red_N_rot']} |")
    return "\n".join(lines)


# ------------------------------------------------------------------ parts

def run_perms(max_n, log):
    out = {}
    for n in range(4, max_n + 1):
        t0 = time.time()
        agg = Agg(n)
        for pi in itertools.permutations(range(n)):
            agg.add(analyze(pi))
        # rebuild worst examples with full per-c tables
        worst_ex = {}
        for k in ("Q_min_minus_d", "loss_residual", "loss_local_at_argmin"):
            worst_ex[k] = [analyze(tuple(p)) for _, p in agg.worst[k][:3]]
        s = agg.summary()
        s["seconds"] = round(time.time() - t0, 1)
        s["coverage"] = f"all {agg.count} permutations of S_{n}, all {n} shifts, both routes, words executed; d from certified table"
        s["worst_tables_md"] = {k: [rule15_table(r) for r in v] for k, v in worst_ex.items()}
        out[n] = s
        log(f"n={n}: {agg.count} perms in {s['seconds']} s; mean losses shift {s['loss_shift']['mean']}, route {s['loss_route_chain']['mean']}, "
            f"junction {s['loss_junction_chain']['mean']}, residual {s['loss_residual']['mean']}, local@argmin {s['loss_local_at_argmin']['mean']}; "
            f"max Q_min-d {s['Q_min_minus_d']['max']}, hist Q_R-B {s['hist_Q_R_minus_B_n']}")
    return out


def reflections(n):
    return [tuple((h - i) % n for i in range(n)) for h in range(n)]


def family_inputs(n, rng):
    fam = {}
    for h in range(n):
        fam[f"reflection_h{h}"] = tuple((h - i) % n for i in range(n))
    for k in range(1, n):
        fam[f"L^{k}_sigma"] = apply_word(sigma(n), "L" * k)
        rot = tuple(sigma(n)[(i + k) % n] for i in range(n))       # sigma o rho_k
        fam[f"sigma_rho{k}"] = rot
    # all single transpositions of sigma_n (small perturbations of a reflection)
    for i in range(n):
        for j in range(i + 1, n):
            p = list(sigma(n)); p[i], p[j] = p[j], p[i]
            fam[f"sigma_swap_{i}_{j}"] = tuple(p)
    if n % 8 == 0:
        b = n // 8
        tau = [0] * n
        for k in range(8):
            for j in range(b):
                tau[b * k + j] = b * ((3 * k) % 8) + j
        fam["C23_family"] = tuple((-tau[i]) % n for i in range(n))
    return fam


def run_families(max_n, log):
    out = {}
    rng = random.Random(0)
    for n in range(4, max_n + 1):
        t0 = time.time()
        fam = family_inputs(n, rng)
        seen = {}
        agg = Agg(n)
        rows = {}
        for name, pi in fam.items():
            if pi in seen:
                rows[name] = {"same_as": seen[pi]}
                continue
            seen[pi] = name
            r = analyze(pi)
            agg.add(r)
            rows[name] = {k: r[k] for k in ("pi", "d", "Q_avg", "Q_min", "argmin_c", "Q_J", "Q_R", "Q_RJ",
                                             "N_X_at_argmin", "N_rot_at_argmin", "loss_shift", "loss_route_chain",
                                             "loss_route_at_argmin", "loss_junction_chain", "loss_local_at_argmin",
                                             "loss_local_min_c", "min_c_local_d_plus_R")}
            rows[name]["loss_residual"] = r.get("loss_residual")
            rows[name]["Q_min_minus_d"] = r.get("Q_min_minus_d")
            rows[name]["proof"] = {k: r["proof"][k] for k in ("s", "x", "T_anchor", "T_min", "g1_s_minus_73",
                                                             "g2_73_minus_max81_82", "g3_max_minus_A_r_n",
                                                             "slack_U_minus_expr34")}
            if name.startswith("reflection") or name == "C23_family":
                rows[name]["rule15_md"] = rule15_table(r)
        s = agg.summary()
        s["seconds"] = round(time.time() - t0, 1)
        s["coverage"] = (f"{agg.count} distinct inputs: all {n} reflections, rotations L^k sigma_n and sigma_n o rho_k, "
                         f"all {n * (n - 1) // 2} single transpositions of sigma_n" + (", C23 family" if n % 8 == 0 else "")
                         + ("; d from certified table" if table(n) else "; NO TABLE: d unknown, residual undefined"))
        s["rows"] = rows
        out[n] = s
        refl = [rows[f"reflection_h{h}"] for h in range(n)]
        log(f"n={n}: {agg.count} inputs in {s['seconds']} s; reflections Q_min-B: {[r['Q_min'] - known.target_diameter(n) for r in refl]}, "
            f"Q_R-B: {[r['Q_R'] - known.target_diameter(n) for r in refl]}, Q_RJ-B: {[r['Q_RJ'] - known.target_diameter(n) for r in refl]}; "
            f"shift loss on sigma_n {rows['reflection_h1']['loss_shift']}, route {rows['reflection_h1']['loss_route_at_argmin']}")
    return out


def run_samples(spec, seed, log):
    out = {}
    for n, k_random, k_perturbed in spec:
        rng = random.Random(seed + n)
        t0 = time.time()
        inputs = {"sigma_n": sigma(n), "rev_n": rev(n)}
        for r in range(k_random):
            q = list(range(n)); rng.shuffle(q); inputs[f"random_{r}"] = tuple(q)
        for r in range(k_perturbed):
            h = rng.randrange(n)
            p = [(h - i) % n for i in range(n)]
            for _ in range(rng.randrange(1, 4)):
                i, j = rng.randrange(n), rng.randrange(n); p[i], p[j] = p[j], p[i]
            inputs[f"perturbed_reflection_{r}"] = tuple(p)
        agg = Agg(n)
        rows = {}
        for name, pi in inputs.items():
            r = analyze(pi)
            agg.add(r)
            rows[name] = {k: r[k] for k in ("pi", "Q_avg", "Q_min", "argmin_c", "Q_J", "Q_R", "Q_RJ",
                                             "N_X_at_argmin", "N_rot_at_argmin", "loss_shift", "loss_route_chain",
                                             "loss_route_at_argmin", "loss_junction_chain")}
            rows[name]["proof"] = {k: r["proof"][k] for k in ("s", "x", "T_anchor", "T_min", "g1_s_minus_73",
                                                             "g2_73_minus_max81_82", "g3_max_minus_A_r_n",
                                                             "slack_U_minus_expr34")}
        s = agg.summary()
        s["seconds"] = round(time.time() - t0, 1)
        s["seed"] = seed + n
        s["coverage"] = f"SAMPLED: {agg.count} inputs (sigma_n, rev_n, {k_random} random, {k_perturbed} perturbed reflections); no table, d unknown; reference B_n"
        s["rows"] = rows
        out[n] = s
        log(f"n={n}: {agg.count} inputs in {s['seconds']} s; mean shift {s['loss_shift']['mean']}, route {s['loss_route_chain']['mean']}, "
            f"junction {s['loss_junction_chain']['mean']}; max Q_min-B {max(int(k) for k in s['hist_Q_min_minus_B_n'])}, "
            f"max Q_R-B {max(int(k) for k in s['hist_Q_R_minus_B_n'])}, slack_U mean {s['slack_U_minus_expr34']['mean']}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perms", type=int, default=0, help="all permutations 4<=n<=PERMS")
    ap.add_argument("--families", type=int, default=0, help="families 4<=n<=FAMILIES")
    ap.add_argument("--samples", action="store_true")
    ap.add_argument("--seed", type=int, default=20260912)
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)
    t0 = time.time()
    logs = []

    def log(msg):
        print(msg, flush=True)
        logs.append(msg)
    meta = {"experiment": EXPERIMENT_VERSION, "construction": su.CONSTRUCTION_VERSION, "core": CORE_VERSION,
            "command": " ".join(sys.argv), "args": vars(args), "python": sys.version.split()[0]}
    if args.perms:
        res = run_perms(args.perms, log)
        with open(os.path.join(OUT_DIR, f"perms{args.tag}.json"), "w") as f:
            json.dump({"meta": meta, "results": res, "log": logs}, f, indent=1, default=to_json)
    if args.families:
        res = run_families(args.families, log)
        with open(os.path.join(OUT_DIR, f"families{args.tag}.json"), "w") as f:
            json.dump({"meta": meta, "results": res, "log": logs}, f, indent=1, default=to_json)
    if args.samples:
        res = run_samples([(20, 30, 10), (50, 20, 10), (100, 10, 6)], args.seed, log)
        with open(os.path.join(OUT_DIR, f"samples{args.tag}.json"), "w") as f:
            json.dump({"meta": meta, "results": res, "log": logs}, f, indent=1, default=to_json)
    log(f"done in {round(time.time() - t0, 1)} s")


if __name__ == "__main__":
    main()
