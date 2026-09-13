"""Extended check of H10 (C27) restricted to the reflection family
pi(i) = h - i mod n, for all h, all c, using the carrier-route variant
(constructions/strict_upper.py, route="carrier") and free reduction.

This extends the exhaustive coverage of C27v beyond n <= 12 (families.json)
for this one structured family, at the cost this module measures. It does
NOT check H10 on general permutations (that needs the full loss_map.py, rule
9 AGENTS.md: n <= 9 exhaustive/n <= 10 with budget, larger n only SAMPLED).

For each n, for each h, computes:
  Q_R  = min_c [2F_c - S_c + R_c]              (route="carrier", raw length)
  Q_RJ = min_c (freely_reduce length at that c)  (min over c, not of the Q_R word)
and reports max_h (Q_R - B_n), max_h (Q_RJ - B_n).

Usage: python3 experiments/reflections_h10_check.py [--max-n 70]
Output: data/runs/reflections_fs/h10_report.json, h10_report.md
"""
import argparse
import json
import os
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "constructions"))
sys.path.insert(0, os.path.join(ROOT, "oracle"))

from strict_upper import shift_word, CONSTRUCTION_VERSION  # noqa: E402
from moves import freely_reduce, CORE_VERSION  # noqa: E402


def bn(n):
    return n * (n - 1) // 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=70)
    ap.add_argument("--min-n", type=int, default=4)
    args = ap.parse_args()

    t0 = time.time()
    per_n = {}
    for n in range(args.min_n, args.max_n + 1):
        tn = time.time()
        B = bn(n)
        worst_raw = None
        worst_red = None
        for h in range(n):
            pi = tuple((h - i) % n for i in range(n))
            Q_R = Q_RJ = None
            for c in range(n):
                sw = shift_word(pi, c, route="carrier")
                L = sw["stats"]["len"]
                red = len(freely_reduce(sw["word"]))
                if Q_R is None or L < Q_R:
                    Q_R = L
                if Q_RJ is None or red < Q_RJ:
                    Q_RJ = red
            if worst_raw is None or (Q_R - B) > worst_raw[0]:
                worst_raw = (Q_R - B, h)
            if worst_red is None or (Q_RJ - B) > worst_red[0]:
                worst_red = (Q_RJ - B, h)
        per_n[n] = {
            "B_n": B,
            "max_QR_minus_Bn": worst_raw[0], "argmax_h_raw": worst_raw[1],
            "max_QRJ_minus_Bn": worst_red[0], "argmax_h_red": worst_red[1],
            "n_mod_4": n % 4,
            "seconds": round(time.time() - tn, 3),
        }
    elapsed = time.time() - t0

    os.makedirs(os.path.join(ROOT, "data", "runs", "reflections_fs"), exist_ok=True)
    summary = {
        "goal": "extend C27v (H10 on reflections) beyond n<=12, all h, all c, exhaustive",
        "core_version": CORE_VERSION, "construction_version": CONSTRUCTION_VERSION,
        "range": f"{args.min_n}<=n<={args.max_n}, all h, all c (exhaustive)",
        "seconds_total": round(elapsed, 2),
        "max_over_all_n_QR_minus_Bn": max(v["max_QR_minus_Bn"] for v in per_n.values()),
        "max_over_all_n_QRJ_minus_Bn": max(v["max_QRJ_minus_Bn"] for v in per_n.values()),
    }
    with open(os.path.join(ROOT, "data", "runs", "reflections_fs", "h10_report.json"), "w") as f:
        json.dump({"summary": summary, "per_n": per_n}, f, indent=1)
    with open(os.path.join(ROOT, "data", "runs", "reflections_fs", "h10_report.md"), "w") as f:
        f.write("# H10 (C27) на отражениях — расширенный диапазон\n\n")
        f.write(f"Диапазон: {summary['range']}. Время: {summary['seconds_total']} с. "
                f"{CORE_VERSION}, {CONSTRUCTION_VERSION}.\n\n")
        f.write(f"max по всем n Q_R - B_n = {summary['max_over_all_n_QR_minus_Bn']}; "
                f"max по всем n Q_RJ - B_n = {summary['max_over_all_n_QRJ_minus_Bn']}.\n\n")
        f.write("| n | n mod 4 | B_n | Q_R-B_n (max по h) | h | Q_RJ-B_n (max по h) | h | сек |\n")
        f.write("|---|---|---|---|---|---|---|---|\n")
        for n in sorted(per_n):
            v = per_n[n]
            f.write(f"| {n} | {v['n_mod_4']} | {v['B_n']} | {v['max_QR_minus_Bn']} | "
                     f"{v['argmax_h_raw']} | {v['max_QRJ_minus_Bn']} | {v['argmax_h_red']} | "
                     f"{v['seconds']} |\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
