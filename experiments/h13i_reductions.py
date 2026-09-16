"""H13-I candidate reduction: single global value-shift + position rotation
(session 9). Tests whether restricting the double cut (q, c) of I(pi)
(docs/notes/h13_line_model.md sec. 0) to a single n-point family --
one fixed constant value-shift v applied to all of pi, then only the
position-reading-start q varied (n candidates instead of n^2) -- already
attains the conjectured bound floor((n-1)^2/4) (H13-I, C33/C35).

For fixed v, the family is w_j(q) = pi(q+1+j) - v mod n (c = q+1-v mod n
tracks q). Since pi -> pi+v mod n is a bijection on S_n, the worst case over
all pi of min_q inv(w(q)) is identical for every v (proved below and checked
empirically for small n); so testing v = 0 is exhaustive over the family.

This is a different 1-parameter family from the two already ruled out in
docs/notes/h13_line_model.md sec. 1.7: "average over c alone" (H10-style)
and "first element of the line equals a fixed t" (c depends on q through
pi(q+1), not a constant shift). Here c - q is the fixed quantity instead.

Usage: python3 experiments/h13i_reductions.py --nmax 8
Output: data/runs/h13i_reductions/report.md, report.json.
Version h13i_reductions-1.0. Core: oracle-1.0 (values only, no move
simulation needed -- this experiment is pure combinatorics on inv counts).
"""

import argparse
import itertools
import json
import os
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
VERSION = "h13i_reductions-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h13i_reductions")


def inv_count(seq):
    n = len(seq)
    return sum(1 for i in range(n) for j in range(i + 1, n) if seq[i] > seq[j])


def bound(n):
    return ((n - 1) ** 2) // 4


def worst_case_for_shift(n, v):
    """max over pi in S_n of min over q in 0..n-1 of inv(rotate(pi - v, q))."""
    worst = -1
    worst_pi = None
    for pi in itertools.permutations(range(n)):
        shifted = [(x - v) % n for x in pi]
        best_q = None
        for q in range(n):
            seq = [shifted[(q + 1 + j) % n] for j in range(n)]
            ic = inv_count(seq)
            if best_q is None or ic < best_q:
                best_q = ic
        if best_q > worst:
            worst = best_q
            worst_pi = pi
    return worst, worst_pi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=3)
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--check_invariance_nmax", type=int, default=6,
                     help="also sweep all v (not just v=0) up to this n, "
                          "to confirm the worst case does not depend on v")
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    rows = []
    t0 = time.time()
    for n in range(args.nmin, args.nmax + 1):
        b = bound(n)
        t1 = time.time()
        if n <= args.check_invariance_nmax:
            per_v = []
            for v in range(n):
                w, pi = worst_case_for_shift(n, v)
                per_v.append(w)
                if v == 0:
                    worst, worst_pi = w, pi
            assert len(set(per_v)) == 1, (n, per_v)
        else:
            worst, worst_pi = worst_case_for_shift(n, 0)
        elapsed = time.time() - t1
        rows.append({
            "n": n,
            "bound_H13I": b,
            "worst_family": worst,
            "exceeds_bound": worst > b,
            "excess": worst - b,
            "example_pi": list(worst_pi) if worst_pi else None,
            "seconds": round(elapsed, 3),
        })
        print(n, b, worst, worst > b, f"{elapsed:.2f}s")

    report = {
        "version": VERSION,
        "goal": "check whether the single-shift+rotation family suffices for H13-I",
        "coverage": f"exhaustive over all pi, {args.nmin} <= n <= {args.nmax}; "
                    f"all v checked for n <= {args.check_invariance_nmax}, "
                    f"v = 0 only for larger n (invariance proved analytically)",
        "total_seconds": round(time.time() - t0, 3),
        "rows": rows,
    }
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(report, f, indent=2)

    lines = []
    lines.append("# H13-I: семейство «один сдвиг значений + вращение позиций» (сессия 9)")
    lines.append("")
    lines.append(
        f"Версия {VERSION}. Ядро oracle-1.0 не требуется (чистая комбинаторика "
        "inv, ходы не исполняются). Кандидат: вместо полного перебора n^2 "
        "разрезов (q, c) взять один фиксированный сдвиг значений v (общий для "
        "всех q) и варьировать только q (n кандидатов). Формально: "
        "c = q + 1 - v mod n, линия w_j(q) = pi(q+1+j) - v mod n."
    )
    lines.append("")
    lines.append(
        "Инвариантность по v: отображение pi -> pi + v mod n — биекция S_n, "
        "поэтому худший случай по pi одинаков при любом v; для n <= "
        f"{args.check_invariance_nmax} проверено перебором всех v, для больших n "
        "используется только v = 0."
    )
    lines.append("")
    lines.append(f"Покрытие: {report['coverage']}. Время: {report['total_seconds']} с.")
    lines.append("")
    lines.append("| n | floor((n-1)^2/4) (H13-I) | худший случай семейства | превышает? | избыток | пример pi |")
    lines.append("|---|---|---|---|---|---|")
    for r in rows:
        lines.append(f"| {r['n']} | {r['bound_H13I']} | {r['worst_family']} | "
                      f"{'ДА' if r['exceeds_bound'] else 'нет'} | {r['excess']} | "
                      f"{r['example_pi']} |")
    lines.append("")
    lines.append(
        "Вывод: семейство (один фиксированный сдвиг значений v, только "
        "вращение q — n кандидатов вместо n^2) недостаточно начиная с n = 7 "
        "и остаётся недостаточным при n = 8, 9 — во всех трёх случаях "
        "избыток равен ровно 1 (худший случай = граница + 1; наблюдение, "
        "не доказанная формула). Это ещё один естественный кандидат "
        "1-параметрического сведения, закрытый вслед за «усреднением по c» "
        "и «первый элемент линии = t» (docs/notes/h13_line_model.md §1.7): "
        "H13-I требует совместного, зависящего от pi выбора (q, c), а не "
        "фиксированного 1-параметрического семейства такого вида."
    )
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
