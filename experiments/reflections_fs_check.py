"""Check of docs/proofs/reflections_FS.md, Lemma 3: exact F_c, S_c on the
reflection family pi(i) = h - i mod n, for all h, all c, as a function of the
parity of m = (h+c) mod n only.

This is a finite check (rule 5/9 AGENTS.md): it does not replace the algebraic
proof in the note, it catches implementation and transcription errors by
comparing the closed formula to constructions/strict_upper.py (strict_upper-1.0)
on every (h, c) pair for n in the given range.

Usage: python3 experiments/reflections_fs_check.py [--max-n 200]
Output: data/runs/reflections_fs/report.md, report.json
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
from moves import CORE_VERSION  # noqa: E402


def P(n):
    return n * n // 4


def predicted_FS(n, m):
    m %= n
    if n % 2 == 1:
        return P(n), 3 * (n - 1) // 2
    if n % 4 == 0:
        return (P(n), 3 * (n // 2 - 1)) if m % 2 == 0 else (P(n), 3 * (n // 2))
    return (P(n) - 1, 3 * (n // 2 - 1)) if m % 2 == 0 else (P(n) + 1, 3 * (n // 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=200)
    ap.add_argument("--min-n", type=int, default=4)
    args = ap.parse_args()

    t0 = time.time()
    mismatches = []
    checked_pairs = 0
    for n in range(args.min_n, args.max_n + 1):
        for h in range(n):
            pi = tuple((h - i) % n for i in range(n))
            for c in range(n):
                sw = shift_word(pi, c, route="n1")
                F_c, S_c = sw["stats"]["F_c"], sw["stats"]["S_c"]
                pf, ps = predicted_FS(n, (h + c) % n)
                checked_pairs += 1
                if (F_c, S_c) != (pf, ps):
                    mismatches.append({"n": n, "h": h, "c": c,
                                        "actual": [F_c, S_c], "predicted": [pf, ps]})
    elapsed = time.time() - t0

    os.makedirs(os.path.join(ROOT, "data", "runs", "reflections_fs"), exist_ok=True)
    report = {
        "goal": "verify docs/proofs/reflections_FS.md Lemma 3 (F_c, S_c on reflections)",
        "core_version": CORE_VERSION, "construction_version": CONSTRUCTION_VERSION,
        "range": f"{args.min_n}<=n<={args.max_n}, all h, all c (exhaustive)",
        "checked_pairs": checked_pairs, "mismatches": len(mismatches),
        "seconds": round(elapsed, 2),
    }
    with open(os.path.join(ROOT, "data", "runs", "reflections_fs", "report.json"), "w") as f:
        json.dump({"summary": report, "mismatches": mismatches[:50]}, f, indent=1)
    with open(os.path.join(ROOT, "data", "runs", "reflections_fs", "report.md"), "w") as f:
        f.write("# Проверка Леммы 3 (docs/proofs/reflections_FS.md)\n\n")
        f.write(f"Диапазон: {report['range']}. Пар (h,c) проверено: {checked_pairs}. "
                f"Расхождений: {len(mismatches)}. Время: {report['seconds']} с. "
                f"{CORE_VERSION}, {CONSTRUCTION_VERSION}.\n")
        if mismatches:
            f.write("\nПервые расхождения:\n")
            for m in mismatches[:20]:
                f.write(f"- n={m['n']} h={m['h']} c={m['c']} "
                        f"actual={m['actual']} predicted={m['predicted']}\n")
    print(json.dumps(report, indent=1))
    if mismatches:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
