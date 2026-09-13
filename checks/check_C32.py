"""Checker for C32: the two-arc reversal word (constructions/reflection_word.py)
sorts every reflection pi_h(i) = h - i mod n in B_n - delta_n(h, 1) moves with
N_X = floor((n-1)^2/4) and N_rot = floor(n^2/4) - delta_n(h, 1).

Parts:
  1. formulas: all h, 4 <= n <= NMAX; each word re-executed with the reference
     moves (oracle-1.0), free reduction checked to be a no-op (no XX, LR, RL);
  2. tables: for 4 <= n <= 10 (certified tables C2v) d(pi_h) == B_n - delta(h, 1)
     for all h (C31 upper half attained; extends C31's range to n = 4, 5);
  3. policy: the horizon policy (candidate_sweep sweep-1.2, mode "horizon",
     H = 3, w = 0.5, thr = 1, c = 0, dir = +1) on pi_{h'} makes exactly the swap
     sequence of policy_like(n, twoA) — both arcs centre-out — for every h' where
     it does not dead-end, 4 <= n <= PMAX; its length then equals
     2 N_X - 2 + (walk between the arcs) + (final walk to c).
Usage: python3 checks/check_C32.py [--nmax 100] [--pmax 24]
Output: data/runs/check_C32/report.json, report.md.  Finite check; the proof
(docs/proofs/C32_two_arc_reversal.md) does not depend on it.
"""

import argparse
import json
import os
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
sys.path.insert(0, os.path.join(ROOT, "exact"))
sys.path.insert(0, os.path.join(ROOT, "constructions"))

from moves import apply_word, freely_reduce, identity, delta, CORE_VERSION  # noqa: E402
from bfs import factorials, rank_perm  # noqa: E402
import reflection_word as rw  # noqa: E402
import candidate_sweep as cs  # noqa: E402

VERSION = "check_C32-1.0"
OUT = os.path.join(ROOT, "data", "runs", "check_C32")


def swap_edges(word, n):
    head, out = 0, []
    for m in word:
        if m == "X":
            out.append(head)
        elif m == "L":
            head = (head + 1) % n
        else:
            head = (head - 1) % n
    return out, head


def part_formulas(nmax):
    rows, fails = [], []
    t0 = time.time()
    for n in range(4, nmax + 1):
        B = n * (n - 1) // 2
        nx_expected = (n - 1) ** 2 // 4
        cases = {}
        for h in range(n):
            r = rw.reflection_word(n, h)
            pi = tuple((h - i) % n for i in range(n))
            w = r["word"]
            ok = (apply_word(pi, w) == identity(n)
                  and freely_reduce(w) == w
                  and r["len"] == B - delta(h, 1, n)
                  and r["N_X"] == nx_expected
                  and r["N_rot"] == n * n // 4 - delta(h, 1, n)
                  and r["len"] == r["N_X"] + r["N_rot"])
            if not ok:
                fails.append({"n": n, "h": h, "len": r["len"], "N_X": r["N_X"], "N_rot": r["N_rot"]})
            cases[r["case"]] = cases.get(r["case"], 0) + 1
        rows.append({"n": n, "B_n": B, "N_X": nx_expected, "cases": cases,
                     "len_sigma": rw.reflection_word(n, 1)["len"],
                     "len_rev": rw.reflection_word(n, n - 1)["len"]})
    return {"nmax": nmax, "seconds": round(time.time() - t0, 1), "fails": fails, "rows": rows}


def part_tables(nmax):
    rows, fails = [], []
    for n in range(4, nmax + 1):
        path = os.path.join(ROOT, "data", "tables", f"dist_n{n}.bin")
        if not os.path.exists(path):
            continue
        with open(path, "rb") as f:
            tab = f.read()
        fact = factorials(n)
        B = n * (n - 1) // 2
        ds = []
        for h in range(n):
            pi = tuple((h - i) % n for i in range(n))
            d = tab[rank_perm(pi, fact)]
            ds.append(d)
            if d != B - delta(h, 1, n):
                fails.append({"n": n, "h": h, "d": d, "formula": B - delta(h, 1, n)})
        rows.append({"n": n, "d_by_h": ds})
    return {"fails": fails, "rows": rows}


def _arc_positions(n, lo, a):
    return {(lo + i) % n for i in range(a)}


def classify_policy(n, twoA, edges):
    """Compare the policy's swap-edge sequence with the two-arc zigzag structure
    about the axis point twoA (first arc centred there, balanced lengths).
    Returns (kind, detail): kind in {"sequential", "interleaved", "other"}."""
    try:
        _, lo1, a, lo2, b = rw.two_arc_edges(n, twoA, "out", "out")
    except ValueError:
        return "other", None
    modes = ("out", "out_mirror", "in", "in_mirror")
    z1 = {m: [e % n for e in rw.zigzag(lo1, a, m)] for m in modes}
    z2 = {m: [e % n for e in rw.zigzag(lo2, b, m)] for m in modes}
    for m1 in modes:
        for m2 in modes:
            if edges == z1[m1] + z2[m2]:
                return "sequential", (m1, m2)
    # interleaved: split edges by arc (an edge e belongs to arc 1 iff e and e+1 are
    # both in arc 1); boundary edges must never be swapped
    p1 = _arc_positions(n, lo1, a)
    s1, s2 = [], []
    for e in edges:
        if e in p1 and (e + 1) % n in p1:
            s1.append(e)
        elif e not in p1 and (e + 1) % n not in p1:
            s2.append(e)
        else:
            return "other", None
    for m1 in modes:
        for m2 in modes:
            if s1 == z1[m1] and s2 == z2[m2]:
                return "interleaved", (m1, m2)
    return "other", None


def part_policy(pmax):
    rows, other, dead = [], [], []
    t0 = time.time()
    for n in range(4, pmax + 1):
        B = n * (n - 1) // 2
        nx = (n - 1) ** 2 // 4
        for hp in range(n):
            pi = tuple((hp - i) % n for i in range(n))
            try:
                r = cs.sweep_word(pi, 0, 1, 1, mode="horizon", horizon=3, weight=0.5)
            except cs.ConstructionError:
                dead.append((n, hp))
                continue
            edges, _ = swap_edges(r["word"], n)
            kind, detail, twoA = "other", None, None
            for cand in (hp, hp + n):
                kind, detail = classify_policy(n, cand, edges)
                if kind != "other":
                    twoA = cand
                    break
            approach = delta(0, edges[0], n)
            final = delta(edges[-1], 0, n)
            junction = sum(delta(edges[i], edges[i + 1], n) - 1 for i in range(len(edges) - 1))
            expl = (r["N_X"] == nx and r["N_rot"] == (r["N_X"] - 1) + approach + final + junction)
            row = {"n": n, "h'": hp, "kind": kind, "modes": detail, "twoA": twoA, "len": r["len"],
                   "len_minus_B": r["len"] - B, "N_X": r["N_X"], "N_rot": r["N_rot"],
                   "approach": approach, "final": final, "junction_extra": junction, "explained": expl}
            rows.append(row)
            if kind == "other":
                other.append({"n": n, "h'": hp, "policy_edges": edges})
    kinds = {}
    for x in rows:
        kinds[x["kind"]] = kinds.get(x["kind"], 0) + 1
    bad = [x for x in rows if not x["explained"]]
    return {"pmax": pmax, "seconds": round(time.time() - t0, 1), "pairs": len(rows), "kinds": kinds,
            "dead_pairs": len(dead), "other": other, "not_explained": bad, "rows": rows}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=100)
    ap.add_argument("--pmax", type=int, default=24)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    f = part_formulas(args.nmax)
    t = part_tables(10)
    p = part_policy(args.pmax)
    ok = not f["fails"] and not t["fails"] and not p["other"] and not p["not_explained"]
    rep = {"version": VERSION, "construction": rw.CONSTRUCTION_VERSION, "candidate": cs.CANDIDATE_VERSION,
           "core": CORE_VERSION, "args": vars(args), "PASS": ok, "formulas": f, "tables": t, "policy": p}
    with open(os.path.join(OUT, "report.json"), "w") as fh:
        json.dump(rep, fh, indent=1)
    lines = [f"# check_C32 — {VERSION}, {rw.CONSTRUCTION_VERSION}, {cs.CANDIDATE_VERSION}, {CORE_VERSION}", "",
             f"Результат: {'PASS' if ok else 'FAIL'}.", "",
             f"1. Формулы: все h при 4 <= n <= {args.nmax} ({sum(r['n'] for r in f['rows'])} слов, {f['seconds']} с): "
             f"слово сортирует (эталон), свободно несократимо, `len = B_n − delta(h,1)`, "
             f"`N_X = floor((n−1)^2/4)`, `N_rot = floor(n^2/4) − delta(h,1)`; нарушений: {len(f['fails'])}.",
             f"2. Таблицы C2v при 4 <= n <= 10: `d(pi_h) = B_n − delta(h,1)` для всех h; нарушений: {len(t['fails'])}.",
             f"3. Политика «горизонт» (c = 0, dir = +1, thr = 1, H = 3, w = 0.5) при 4 <= n <= {args.pmax}: "
             f"{p['pairs']} пар (n, h') без тупика ({p['dead_pairs']} тупиков): обмены — два зигзага дуг "
             f"{p['kinds']} (sequential: конкатенация зигзагов дуг; interleaved: зигзаги дуг чередуются, "
             f"граничные рёбра не трогаются; other: иное — {len(p['other'])}); "
             f"`N_X = floor((n−1)^2/4)` и `N_rot = N_X − 1 + подход + финал + лишние шаги на стыках` "
             f"во всех, кроме {len(p['not_explained'])}; {p['seconds']} с.", "",
             "| n | B_n | N_X | len(sigma_n) | len(rev_n) | случаи |", "|---|---|---|---|---|---|"]
    for r in f["rows"]:
        if r["n"] <= 16 or r["n"] % 10 == 0:
            lines.append(f"| {r['n']} | {r['B_n']} | {r['N_X']} | {r['len_sigma']} | {r['len_rev']} | {r['cases']} |")
    lines += ["", "Политика по h' = 0..n−1: `len − B_n` (вид; подход + финал + стыки); `—` тупик:", ""]
    for n in range(4, args.pmax + 1):
        by = {x["h'"]: x for x in p["rows"] if x["n"] == n}
        cells = []
        for hp in range(n):
            x = by.get(hp)
            cells.append("—" if x is None else
                         f"{x['len_minus_B']:+d} ({x['kind'][0]}; {x['approach']}+{x['final']}+{x['junction_extra']})")
        lines.append(f"- n = {n}: " + ", ".join(cells))
    with open(os.path.join(OUT, "report.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines[:8]))
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
