"""Structure of geodesics (shortest sorting words) on hard inputs, n <= 10.

Question (PLAN section 8, session 5): construction N1 sorts the cycles of
f_c(i) = pi(i) + c one at a time; on affine inputs pi(i) = a*i + b this loses
Theta(n) against B_n (H11).  What do the geodesics do instead?

For a sorting word w of pi (executed left to right, head starts at circle
position 0) the final head position c(w) = (#L - #R) mod n is the shift of
the word: element pi[i] (initially at circle position i) ends at circle
position pi[i] + c.  So every geodesic has its own shift, and its swaps can be
classified against the cycles of f_c for that c:

  * an X at circle edge {h, h+1} moves the element at h by +1 and the element
    at h+1 by -1.  A move is USEFUL for an element if it decreases its
    remaining circular distance to its target (ties at distance n/2 count as
    useful).  A swap is double (both useful), single, or bad.
    Potential Phi = sum of remaining distances: double -2, single 0, bad +2,
    hence N_X = F_c/2 + #single + 2*#bad exactly (F_c = Phi at the start).
  * a double swap is CROSS if its two elements start in different cycles of
    f_c (or one of them is a fixed point of f_c).  Construction N1 plans
    only pairs inside one cycle (2F(C) - S(C) letters per cycle); a cross
    double swap can still happen there by accident, when a bystander of
    another cycle is pushed towards its own target.  The count is reported
    for the construction's words as well.

Edge-load lower bound (pairing bound), per shift c: for every element take
its shortest arc to the target (antipodal arc positive, rule 16); p_v, m_v =
number of arcs crossing edge v in the + / - direction.  Every X on edge v
serves at most one + crossing and one - crossing, so
    N_X >= LBX(c) = sum_v max(p_v, m_v)
for any word with shift c whose elements travel along shortest arcs; the
per-word count is exact when no element deviates.  Reported next to the
construction's N_X and the geodesics' N_X.

Geodesic DAG from the certified table dist_n<n>.bin: a state s reached from
pi in t moves lies on a geodesic iff d(s) = d(pi) - t.  Counting is exact
(big ints), sampling is uniform over all geodesics.

Usage: python3 experiments/geodesic_structure.py --n 10 [--samples 300]
       [--inputs affine|all] [--seed 1]
Output: data/runs/geodesic_structure/report_n<n>.json and .md
Version geodesic_structure-1.0; core oracle-1.0; construction strict_upper-1.0.
"""

import argparse
import json
import os
import random
import sys
import time
from collections import Counter
from math import gcd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
sys.path.insert(0, os.path.join(ROOT, "exact"))
sys.path.insert(0, os.path.join(ROOT, "constructions"))

from moves import apply_move, apply_word, identity, sigma, delta, CORE_VERSION  # noqa: E402
from bfs import factorials, rank_perm  # noqa: E402
import strict_upper as su  # noqa: E402

VERSION = "geodesic_structure-1.0"
OUT_DIR = os.path.join(ROOT, "data", "runs", "geodesic_structure")


# ----------------------------------------------------------------- table

def load_table(n):
    path = os.path.join(ROOT, "data", "tables", f"dist_n{n}.bin")
    with open(path, "rb") as f:
        return f.read()


class Dist:
    def __init__(self, n):
        self.n = n
        self.fact = factorials(n)
        self.tab = load_table(n)
        assert len(self.tab) == self.fact[n]

    def __call__(self, p):
        return self.tab[rank_perm(p, self.fact)]


# ----------------------------------------------------------------- geodesic DAG

def geodesic_dag(pi, dist):
    """Layers of the geodesic DAG pi -> id: list of dicts state -> (#paths from pi).
    Also returns, per state, the number of completions to id (for uniform sampling)."""
    n = len(pi)
    d0 = dist(pi)
    layers = [{pi: 1}]
    for t in range(1, d0 + 1):
        nxt = {}
        want = d0 - t
        for s, cnt in layers[-1].items():
            for m in ("L", "R", "X"):
                u = apply_move(s, m)
                if dist(u) == want:
                    nxt[u] = nxt.get(u, 0) + cnt
        layers.append(nxt)
    assert list(layers[-1]) == [identity(n)]
    # completions
    comp = [dict() for _ in layers]
    comp[-1][identity(n)] = 1
    for t in range(d0 - 1, -1, -1):
        for s in layers[t]:
            tot = 0
            for m in ("L", "R", "X"):
                u = apply_move(s, m)
                if u in comp[t + 1]:
                    tot += comp[t + 1][u]
            comp[t][s] = tot
    return layers, comp


def sample_geodesic(pi, layers, comp, rng):
    """Uniform random geodesic word from pi to id."""
    s, word = pi, []
    for t in range(len(layers) - 1):
        opts = []
        for m in ("L", "R", "X"):
            u = apply_move(s, m)
            if u in comp[t + 1]:
                opts.append((m, u, comp[t + 1][u]))
        tot = sum(o[2] for o in opts)
        r = rng.randrange(tot)
        for m, u, c in opts:
            if r < c:
                word.append(m)
                s = u
                break
            r -= c
    return "".join(word)


def nx_distribution(pi, layers):
    """Exact distribution of N_X over all geodesics (DP with polynomial in x)."""
    n = len(pi)
    poly = {pi: {0: 1}}
    for t in range(len(layers) - 1):
        nxt = {}
        for s, pl in poly.items():
            for m in ("L", "R", "X"):
                u = apply_move(s, m)
                if u in layers[t + 1]:
                    acc = nxt.setdefault(u, {})
                    for k, v in pl.items():
                        kk = k + (1 if m == "X" else 0)
                        acc[kk] = acc.get(kk, 0) + v
        poly = nxt
    return dict(sorted(poly[identity(n)].items()))


# ----------------------------------------------------------------- swap classification

def shift_of(word, n):
    return (word.count("L") - word.count("R")) % n


def elements_by_cycle(pi, c):
    """Map element -> cycle id of f_c (fixed points get their own negative ids)."""
    n = len(pi)
    f = [(pi[i] + c) % n for i in range(n)]
    cyc_id = {}
    for k, cyc in enumerate(su.cycles_of(f)):
        for pos in cyc:
            cyc_id[pi[pos]] = k
    for i in range(n):
        if f[i] == i:
            cyc_id[pi[i]] = -1 - i
    return cyc_id, f


def edge_loads(pi, c):
    """p_v, m_v: crossings of edge {v, v+1} by shortest arcs (antipodal positive)."""
    n = len(pi)
    p = [0] * n
    m = [0] * n
    for i in range(n):
        d = su.signed_step(i, (pi[i] + c) % n, n)
        for v in su.arc_edges(i, d, n):
            if d > 0:
                p[v] += 1
            else:
                m[v] += 1
    return p, m


def classify_word(pi, word):
    """Walk the word in the circle model; classify every X."""
    n = len(pi)
    c = shift_of(word, n)
    cyc_id, f = elements_by_cycle(pi, c)
    circle = list(pi)                         # circle[pos] = element
    head = 0
    target = {pi[i]: (pi[i] + c) % n for i in range(n)}
    cnt = Counter()
    phi = sum(delta(i, target[pi[i]], n) for i in range(n))
    assert phi == su_F(pi, c)
    for ch in word:
        if ch == "L":
            head = (head + 1) % n
        elif ch == "R":
            head = (head - 1) % n
        else:
            a, b = head, (head + 1) % n
            ea, eb = circle[a], circle[b]
            ua = delta(b, target[ea], n) <= delta(a, target[ea], n) and delta(a, target[ea], n) > 0
            ub = delta(a, target[eb], n) <= delta(b, target[eb], n) and delta(b, target[eb], n) > 0
            # useful = strictly closer, or equal (antipodal tie) but not already home
            k = int(ua) + int(ub)
            kind = {2: "double", 1: "single", 0: "bad"}[k]
            cnt[kind] += 1
            if k == 2:
                same = cyc_id[ea] == cyc_id[eb] and cyc_id[ea] >= 0
                cnt["double_same" if same else "double_cross"] += 1
            circle[a], circle[b] = eb, ea
    assert all(circle[target[e]] == e for e in range(n))
    assert apply_word(tuple(pi), word) == identity(n)
    return {"c": c, "N_X": word.count("X"), "N_rot": len(word) - word.count("X"),
            "double": cnt["double"], "single": cnt["single"], "bad": cnt["bad"],
            "double_same": cnt["double_same"], "double_cross": cnt["double_cross"],
            "n_cycles": len(su.cycles_of(f))}


def su_F(pi, c):
    n = len(pi)
    return sum(delta(i, (pi[i] + c) % n, n) for i in range(n))


# ----------------------------------------------------------------- construction side

def construction_stats(pi):
    """Per shift: N1 carrier-route word statistics + pairing bound."""
    n = len(pi)
    rows = []
    for c in range(n):
        sw = su.shift_word(pi, c, route="carrier")
        st = sw["stats"]
        cl = classify_word(pi, sw["word"])
        p, m = edge_loads(pi, c)
        lbx = sum(max(a, b) for a, b in zip(p, m))
        rows.append({"c": c, "K": st["n_cycles"], "F": st["F_c"], "S": st["S_c"],
                     "R": st["route_len"], "len": st["len"], "N_X": st["N_X"],
                     "N_rot": st["N_rot"], "LBX": lbx, "F_half": st["F_c"] / 2,
                     "double": cl["double"], "single": cl["single"], "bad": cl["bad"],
                     "double_cross": cl["double_cross"]})
    return rows


# ----------------------------------------------------------------- inputs

def affine_inputs(n):
    out = []
    for a in range(2, n):
        if gcd(a, n) != 1:
            continue
        for b in range(n):
            out.append((f"affine a={a} b={b}", tuple((a * i + b) % n for i in range(n))))
    return out


def worst_affine(n, dist):
    """Affine inputs with a != n-1 (reflections are listed separately): the three
    with the largest min_c len_R - B_n and the three with the largest loss
    min_c len_R - d(pi) of the construction against the true distance."""
    B = n * (n - 1) // 2
    scored = []
    for name, pi in affine_inputs(n):
        if name.startswith(f"affine a={n - 1} "):
            continue
        best = min(su.shift_word(pi, c, route="carrier")["stats"]["len"] for c in range(n))
        scored.append((best - B, best - dist(pi), dist(pi), name, pi))
    by_excess = sorted(scored, key=lambda r: (-r[0], -r[1]))[:3]
    by_loss = sorted(scored, key=lambda r: (-r[1], -r[0]))[:3]
    out, seen = [], set()
    for ex, loss, d, nm, pi in by_excess + by_loss:
        if pi in seen:
            continue
        seen.add(pi)
        out.append((nm + f" (min_c len_R - B_n = {ex}, - d = {loss}, d = {d})", pi))
    return out


def named_inputs(n, dist, which):
    inputs = [("sigma_n", sigma(n))]
    if n >= 6:
        h = 2
        inputs.append((f"reflection h={h}", tuple((h - i) % n for i in range(n))))
    inputs += worst_affine(n, dist)
    if n == 10:
        inputs.append(("C27 part-2 counterexample (9,2,5,0,1,8,7,6,3,4)",
                       (9, 2, 5, 0, 1, 8, 7, 6, 3, 4)))
    if which == "all":
        inputs += affine_inputs(n)
    return inputs


# ----------------------------------------------------------------- main

def analyse(pi, dist, samples, rng):
    n = len(pi)
    t0 = time.time()
    layers, comp = geodesic_dag(pi, dist)
    n_geo = comp[0][pi]
    nxd = nx_distribution(pi, layers)
    words = [sample_geodesic(pi, layers, comp, rng) for _ in range(samples)]
    cls = [classify_word(pi, w) for w in words]
    by_c = Counter(x["c"] for x in cls)
    cons = construction_stats(pi)
    best_cons = min(cons, key=lambda r: r["len"])
    # geodesic summary per shift
    per_c = {}
    for x in cls:
        acc = per_c.setdefault(x["c"], Counter())
        acc["count"] += 1
        for k in ("N_X", "N_rot", "double", "single", "bad", "double_same", "double_cross"):
            acc[k] += x[k]
    per_c_rows = []
    for c in sorted(per_c):
        acc = per_c[c]
        k = acc["count"]
        row = {"c": c, "geodesics_sampled": k}
        for key in ("N_X", "N_rot", "double", "single", "bad", "double_same", "double_cross"):
            row[key + "_mean"] = round(acc[key] / k, 2)
        crow = cons[c]
        row.update({"F": crow["F"], "LBX": crow["LBX"], "K": crow["K"],
                    "construction_len": crow["len"], "construction_N_X": crow["N_X"],
                    "construction_double_cross": crow["double_cross"]})
        per_c_rows.append(row)
    # min over sampled geodesics of cross double swaps
    min_cross = min(x["double_cross"] for x in cls)
    max_cross = max(x["double_cross"] for x in cls)
    zero_cross = sum(1 for x in cls if x["double_cross"] == 0)
    example = min(zip(cls, words), key=lambda cw: (cw[0]["double_cross"], cw[0]["N_X"]))
    return {
        "pi": list(pi), "d": dist(pi), "B_n": n * (n - 1) // 2,
        "geodesics": n_geo, "dag_states": sum(len(l) for l in layers),
        "N_X_distribution": {str(k): v for k, v in nxd.items()},
        "geodesic_shifts": {str(k): v for k, v in sorted(by_c.items())},
        "min_double_cross_sampled": min_cross, "max_double_cross_sampled": max_cross,
        "sampled_with_zero_cross": zero_cross,
        "example_min_cross": {"word": example[1], **example[0]},
        "construction_best": best_cons,
        "construction_min_LBX_plus": min(r["LBX"] for r in cons),
        "per_shift": per_c_rows,
        "construction_per_shift": cons,
        "seconds": round(time.time() - t0, 2),
    }


def write_md(n, results, args, path):
    B = n * (n - 1) // 2
    L = [f"# Геодезические на трудных входах, n = {n}", "",
         f"Версии: {VERSION}, {su.CONSTRUCTION_VERSION}, {CORE_VERSION}; таблица "
         f"`data/tables/dist_n{n}.bin` (C2v). Команда: `python3 experiments/geodesic_structure.py "
         f"--n {n} --samples {args.samples} --inputs {args.inputs} --seed {args.seed}`.", "",
         "Обозначения: c — сдвиг слова (конечное положение головы); F — сумма кратчайших дуг при этом c; "
         "LBX — оценка спаривания `sum_v max(p_v, m_v)`; double/single/bad — обмены, полезные для двух / "
         "одного / нуля элементов; cross — двойные обмены между элементами разных циклов f_c "
         "(конструкция N1 их не делает). `N_X = F/2 + single + 2·bad` — тождество.", ""]
    L.append("| вход | d | B_n | геодезич. | N_X (все геодезические) | сдвиги геодезических | "
             "N1: c*, len, N_X, cross | геодезические (среднее по выборке при c*): N_X, double_same, cross, single, bad | LBX при c* | min cross по выборке |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for name, r in results:
        bc = r["construction_best"]
        pc = next((x for x in r["per_shift"] if x["c"] == bc["c"]), None)
        geo = (f"{pc['N_X_mean']}, {pc['double_same_mean']}, {pc['double_cross_mean']}, "
               f"{pc['single_mean']}, {pc['bad_mean']}" if pc else "нет геодезических с этим c")
        L.append(f"| {name} | {r['d']} | {B} | {r['geodesics']} | {r['N_X_distribution']} | "
                 f"{r['geodesic_shifts']} | {bc['c']}, {bc['len']}, {bc['N_X']}, {bc['double_cross']} | "
                 f"{geo} | {bc['LBX']} | {r['min_double_cross_sampled']} "
                 f"(из {args.samples}; с нулём: {r['sampled_with_zero_cross']}) |")
    L.append("")
    for name, r in results:
        L.append(f"## {name}")
        L.append("")
        L.append(f"pi = {tuple(r['pi'])}, d = {r['d']}, геодезических {r['geodesics']}, "
                 f"состояний ДАГа {r['dag_states']}, {r['seconds']} с.")
        L.append("")
        L.append("| c | K | F | LBX | N1 len | N1 N_X | N1 cross | геод. в выборке | геод. N_X | double_same | cross | single | bad | N_rot |")
        L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        geo_by_c = {x["c"]: x for x in r["per_shift"]}
        for cr in r["construction_per_shift"]:
            g = geo_by_c.get(cr["c"])
            gs = (f"{g['geodesics_sampled']} | {g['N_X_mean']} | {g['double_same_mean']} | {g['double_cross_mean']} | "
                  f"{g['single_mean']} | {g['bad_mean']} | {g['N_rot_mean']}") if g else "0 | | | | | | "
            L.append(f"| {cr['c']} | {cr['K']} | {cr['F']} | {cr['LBX']} | {cr['len']} | {cr['N_X']} | "
                     f"{cr['double_cross']} | {gs} |")
        ex = r["example_min_cross"]
        L.append("")
        L.append(f"Пример геодезической с минимумом cross в выборке: c = {ex['c']}, N_X = {ex['N_X']}, "
                 f"cross = {ex['double_cross']}, single = {ex['single']}, bad = {ex['bad']}: `{ex['word']}`")
        L.append("")
    with open(path, "w") as f:
        f.write("\n".join(L))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=10)
    ap.add_argument("--samples", type=int, default=300)
    ap.add_argument("--inputs", default="affine", choices=["affine", "all"])
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    n = args.n
    rng = random.Random(args.seed)
    dist = Dist(n)
    os.makedirs(OUT_DIR, exist_ok=True)
    results = []
    t0 = time.time()
    for name, pi in named_inputs(n, dist, args.inputs):
        r = analyse(pi, dist, args.samples, rng)
        results.append((name, r))
        bc = r["construction_best"]
        print(f"{name}: d={r['d']} geodesics={r['geodesics']} N_X={r['N_X_distribution']} "
              f"shifts={r['geodesic_shifts']} N1 best c={bc['c']} len={bc['len']} "
              f"min_cross={r['min_double_cross_sampled']} ({r['seconds']} s)", flush=True)
    meta = {"version": VERSION, "construction": su.CONSTRUCTION_VERSION, "core": CORE_VERSION,
            "n": n, "samples": args.samples, "seed": args.seed, "inputs": args.inputs,
            "seconds": round(time.time() - t0, 1),
            "command": f"python3 experiments/geodesic_structure.py --n {n} --samples {args.samples} "
                       f"--inputs {args.inputs} --seed {args.seed}"}
    with open(os.path.join(OUT_DIR, f"report_n{n}.json"), "w") as f:
        json.dump({"meta": meta, "results": [{"name": nm, **r} for nm, r in results]}, f,
                  ensure_ascii=False, indent=1)
    write_md(n, results, args, os.path.join(OUT_DIR, f"report_n{n}.md"))
    print("done", meta["seconds"], "s")


if __name__ == "__main__":
    main()
