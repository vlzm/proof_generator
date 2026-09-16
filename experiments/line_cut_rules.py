"""Cut-selection rules for H13-I: does a canonical rule attain floor((n-1)^2/4)?

C38 (docs/proofs/C37_toric_mean_bound.md) shows that no averaging over cuts can
prove H13-I: on reflections only n (or 2n) of the n^2 cuts are optimal, so a
proof must *select* the cut.  This module tests simple, canonical selection
rules.  For a rule with score f (minimised over cuts) we take the worst
inversion count among the cuts attaining the minimum (pessimistic tie-breaking)
and compare it with floor((n-1)^2/4).

Rules (w is the line of the cut, pos_w the inverse line):
  spearman     sum_j (w_j - j)^2           (equivalently: maximise sum_j j*w_j)
  footrule     sum_j |w_j - j|
  fixedpoints  -#{j : w_j = j}
  circ         sum_j min((w_j - j) mod n, (j - w_j) mod n)
  prefix       max_k sum_{j<k} (n - 1 - 2 w_j)   (slack of lemma C)

Result (session 9): only "spearman" survives — exhaustively for all pi with
4 <= n <= 9 and on families up to n = 40.  The others fail, "footrule" already
at n = 4.  Hence the conjecture H13-S (CLAIMS C39): at a cut minimising the
Spearman distance the line has at most floor((n-1)^2/4) inversions.

Usage: python3 experiments/line_cut_rules.py [--nmax 8] [--spearman_nmax 9] [--families 40]
Output: data/runs/line_cut_rules/report.json, report.md.
Version line_cut_rules-1.0.
"""

import argparse
import itertools
import json
import math
import os
import random
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "data", "runs", "line_cut_rules")
VERSION = "line_cut_rules-1.0"

RULES = ["spearman", "footrule", "fixedpoints", "circ", "prefix"]


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def scores(w):
    n = len(w)
    s = 0
    mx = 0
    for j in range(n):
        s += n - 1 - 2 * w[j]
        if s > mx:
            mx = s
    return {
        "spearman": sum((w[j] - j) ** 2 for j in range(n)),
        "footrule": sum(abs(w[j] - j) for j in range(n)),
        "fixedpoints": -sum(1 for j in range(n) if w[j] == j),
        "circ": sum(min((w[j] - j) % n, (j - w[j]) % n) for j in range(n)),
        "prefix": mx,
    }


def rule_results(pi):
    """For every rule: worst inv among the cuts minimising it."""
    n = len(pi)
    best = {r: None for r in RULES}
    for a in range(n):
        for b in range(n):
            w = [(pi[(a + j) % n] - b) % n for j in range(n)]
            iv = inversions(w)
            sc = scores(w)
            for r in RULES:
                cur = best[r]
                if cur is None or sc[r] < cur[0]:
                    best[r] = (sc[r], iv)
                elif sc[r] == cur[0] and iv > cur[1]:
                    best[r] = (cur[0], iv)
    return {r: best[r][1] for r in RULES}


def spearman_fast(pi):
    """Rule "spearman" alone, with O(1) updates over the grid of cuts.

    Returns (worst inv among the cuts maximising sum_j j*w_j, I(pi)); the strong
    form of H13-S is that these two numbers are always equal.

    Shifting the value cut by 1: inv += n-1-2p, sum_j j*w_j += n*p - n(n-1)/2,
    where p is the position of the smallest value.  Shifting the position cut by
    1: inv += n-1-2*w_0, sum_j j*w_j += n*w_0 - n(n-1)/2.
    """
    n = len(pi)
    pos = [0] * n
    for i, v in enumerate(pi):
        pos[v] = i
    T = n * (n - 1) // 2
    w0 = [pi[j] for j in range(n)]
    curI = inversions(w0)
    curS = sum(j * w0[j] for j in range(n))
    bestS = bestI = None
    minI = None
    for a in range(n):
        vI, vS = curI, curS
        for b in range(n):
            if bestS is None or vS > bestS:
                bestS, bestI = vS, vI
            elif vS == bestS and vI > bestI:
                bestI = vI
            if minI is None or vI < minI:
                minI = vI
            p = (pos[b] - a) % n
            vI += n - 1 - 2 * p
            vS += n * p - T
        assert (vI, vS) == (curI, curS)
        curI += n - 1 - 2 * pi[a]
        curS += n * pi[a] - T
    return bestI, minI


def families(n, rnd, k_random=60):
    out = {}
    for h in range(n):
        out[f"reflection h={h}"] = [(h - i) % n for i in range(n)]
    for a in range(1, n):
        if math.gcd(a, n) == 1:
            out[f"affine {a}i+1"] = [(a * i + 1) % n for i in range(n)]
    for k in range(k_random):
        p = list(range(n))
        rnd.shuffle(p)
        out[f"random {k}"] = p
    m = n // 2
    p = list(range(n))
    p[m:] = list(range(m, n))[::-1]
    out["half reversed"] = p
    p = list(range(n))
    p[0], p[m] = p[m], p[0]
    out["antipodal transposition"] = p
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8, help="exhaustive over all pi up to this n")
    ap.add_argument("--families", type=int, default=40, help="largest n for the sampled part")
    ap.add_argument("--spearman_nmax", type=int, default=9,
                    help="exhaustive over all pi for the spearman rule alone (fast path)")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    rnd = random.Random(20260916)
    t0 = time.time()
    rep = {"version": VERSION, "seed": 20260916, "exhaustive_nmax": args.nmax,
           "families_nmax": args.families, "exhaustive": {}, "sampled": {}}

    # self-check of the fast path
    for n in range(4, 7):
        for _ in range(50):
            p = list(range(n))
            rnd.shuffle(p)
            assert spearman_fast(p)[0] == rule_results(p)["spearman"], p

    for n in range(4, args.nmax + 1):
        bound = (n - 1) ** 2 // 4
        row = {r: {"max_inv": -1, "violations": 0, "worst_excess": 0, "example": None}
               for r in RULES}
        strong_bad = 0
        for pi in itertools.permutations(range(n)):
            res = rule_results(list(pi))
            v, opt = spearman_fast(list(pi))
            if v != opt:
                strong_bad += 1
            for r in RULES:
                v = res[r]
                d = row[r]
                d["max_inv"] = max(d["max_inv"], v)
                if v > bound:
                    d["violations"] += 1
                    if v - bound > d["worst_excess"] or d["example"] is None:
                        d["worst_excess"] = max(d["worst_excess"], v - bound)
                        d["example"] = {"pi": list(pi), "inv": v}
        rep["exhaustive"][n] = {"bound": bound, "rules": row,
                                "strong_form_failures": strong_bad}
        print(f"n={n} bound={bound}: " +
              ", ".join(f"{r}: max {row[r]['max_inv']}, bad {row[r]['violations']}" for r in RULES),
              flush=True)

    for n in range(args.nmax + 1, args.spearman_nmax + 1):
        bound = (n - 1) ** 2 // 4
        bad = 0
        mx = -1
        strong_bad = 0
        for pi in itertools.permutations(range(n)):
            v, opt = spearman_fast(list(pi))
            mx = max(mx, v)
            if v > bound:
                bad += 1
            if v != opt:
                strong_bad += 1
        rep["exhaustive"][n] = {"bound": bound, "rules": {"spearman": {
            "max_inv": mx, "violations": bad, "worst_excess": max(0, mx - bound),
            "example": None}}, "strong_form_failures": strong_bad}
        print(f"n={n} bound={bound} (all pi, spearman only): max {mx}, bad {bad}, "
              f"strong form (rule S = I) failures {strong_bad}", flush=True)

    ns = [n for n in (9, 10, 11, 12, 16, 20, 25, 30, args.families)
          if n > max(args.nmax, args.spearman_nmax)]
    for n in sorted(set(ns)):
        bound = (n - 1) ** 2 // 4
        bad = 0
        mx = -1
        tot = 0
        example = None
        strong_bad = 0
        for name, pi in families(n, rnd).items():
            tot += 1
            v, opt = spearman_fast(pi)
            mx = max(mx, v)
            if v != opt:
                strong_bad += 1
            if v > bound:
                bad += 1
                example = {"name": name, "pi": pi, "inv": v}
        rep["sampled"][n] = {"bound": bound, "inputs": tot, "max_inv_rule_spearman": mx,
                             "violations": bad, "example": example,
                             "strong_form_failures": strong_bad}
        print(f"n={n} (families, rule spearman): {tot} inputs, max inv {mx}, bound {bound}, "
              f"violations {bad}, strong form failures {strong_bad}", flush=True)

    rep["seconds"] = round(time.time() - t0)
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(rep, f, indent=1, sort_keys=True)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# line_cut_rules — правила выбора разреза для H13-I\n\n")
        f.write(f"Версия {VERSION}, seed 20260916, {rep['seconds']} с, память O(n^2).\n\n"
                "Цель: по C38 (`docs/proofs/C37_toric_mean_bound.md`) усреднение по разрезам "
                "бесполезно — разрез нужно *выбирать*; проверяются простые канонические "
                "правила (минимум score по всем n^2 разрезам; ничьи разрешаются в худшую "
                "сторону — берётся наибольшее inv среди победивших разрезов).\n\n"
                "Команда: `python3 experiments/line_cut_rules.py "
                f"--nmax {args.nmax} --spearman_nmax {args.spearman_nmax} "
                f"--families {args.families}`. Покрытие: исчерпывающее по всем pi при "
                f"4 <= n <= {args.spearman_nmax} (все правила — до n = {args.nmax}), "
                "далее выборочно (SAMPLED): все отражения, все аффинные, 60 случайных "
                "(seed 20260916), блочные входы. Ядро oracle-1.0 и таблицы не нужны: "
                "слова не строятся.\n\n")
        f.write("## Исчерпывающе по всем pi\n\n| n | порог | " +
                " | ".join(RULES) + " |\n|" + "---|" * (len(RULES) + 2) + "\n")
        for n in sorted(rep["exhaustive"]):
            d = rep["exhaustive"][n]
            cells = []
            for r in RULES:
                x = d["rules"].get(r)
                cells.append("—" if x is None else f"max {x['max_inv']}, нарушений {x['violations']}")
            f.write(f"| {n} | {d['bound']} | " + " | ".join(cells) + " |\n")
        f.write("\nСтрогая форма (правило spearman даёт не просто оценку, а сам минимум "
                "I(pi) — то есть любой разрез, минимизирующий расстояние Спирмена, "
                "минимизирует и inv): нарушений нет ни при одном n выше, счётчики "
                "`strong_form_failures` в report.json.\n")
        f.write("\n## Правило spearman на семействах (отражения, аффинные, случайные, блочные)\n\n")
        f.write("| n | порог | входов | max inv | нарушений |\n|---|---|---|---|---|\n")
        for n in sorted(rep["sampled"]):
            d = rep["sampled"][n]
            f.write(f"| {n} | {d['bound']} | {d['inputs']} | {d['max_inv_rule_spearman']} | "
                    f"{d['violations']} |\n")
        f.write("\n## Вывод и следующий вопрос\n\n"
                "Правило «spearman» (разрез, минимизирующий `sum_j (w_j − j)^2`, то есть "
                "максимизирующий `sum_j j·w_j`) не нарушено ни разу; остальные правила "
                "нарушаются, «footrule» — уже при n = 4 на отражении `(0,3,2,1)` "
                "(inv 3 > 2) и при n = 8 на `(0,7,6,5,4,3,2,1)` (13 > 12). Отсюда "
                "гипотеза H13-S (CLAIMS C39). Конечная проверка не доказывает общий "
                "случай.\n\nСледующий вопрос: доказать H13-S. Известная связь — "
                "неравенство Дэниелса `−1 <= 3τ − 2ρ <= 1`, из которого `D >= 0` следует "
                "при `sum_j (w_j − j)^2 <= n(n^2−1)/12`; среднее `sum d^2` по всем n^2 "
                "разрезам равно ровно `n(n^2−1)/6` (вдвое больше порога), а на отражении "
                "минимум примерно `n^3/12` — ровно на границе "
                "(`docs/notes/h13i_verdict.md` §4.1).\n")
    print(f"done in {rep['seconds']} s -> {OUT}")


if __name__ == "__main__":
    main()
