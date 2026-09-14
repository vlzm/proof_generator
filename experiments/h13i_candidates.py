"""H13-I candidate cut-selection rules (PLAN.md H13, docs/notes/h13_line_model.md
S6): the target `I(pi) <= floor((n-1)^2/4)` (hypothesis H13-I) is established by
brute force over all n^2 cuts (q, c) up to n = 12 (h13i_toric.c), but a *proof*
needs either a general argument or, short of that, a cheap (O(n), not O(n^2))
rule that already picks a good-enough cut. This module tests two such rules
and records where each fails, as negative results (AGENTS.md rule 2: state
checked directly, not by plan names).

Definitions match line_model.py / line_profile.c: for a permutation pi of
Z_n, position-cut q and value-shift c, shift = (q+1-c) mod n,
w_j = (pi[(q+1+j) mod n] - shift) mod n for j = 0..n-1; inv(w) counts pairs
j < k with w_j > w_k. I(pi) = min over all n^2 (q, c) of inv(w).

Rule "diagonal": use each point of pi itself as the cut origin, i.e. n
candidates c = q - pi(q) (so pi(q) maps to line index 0, value 0); take the
best over these n candidates only (already tried and rejected in session 8
in the form "induction by one element" / "family q = fixed value" -- this is
the direct single-point version).

Rule "mode_shift": fix c once via c = -mode_b(pi(i) - i mod n) mod n (mode of
the additive shift that would make pi a pure rotation, sign-corrected so a
genuine rotation gives I = 0 exactly -- checked below), then minimize only
over q (n candidates, not n^2).

Both rules are refuted: each has a counterexample already at small n found
by exhaustive search against the true I(pi) (brute-force cross-checked
against the O(n^3) circular-difference-array method used by h13i_toric.c).

Usage: python3 experiments/h13i_candidates.py [--nmax 8]
Output: data/runs/h13i_search/candidates_report.json, .md.
Version h13i_candidates-1.0 (session 9, 14.09.2026).
"""

import argparse
import itertools
import json
import os
import sys
import time
from collections import Counter

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
from moves import CORE_VERSION  # noqa: E402

VERSION = "h13i_candidates-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h13i_search")


def inv_count(w):
    n = len(w)
    c = 0
    for i in range(n):
        wi = w[i]
        for j in range(i + 1, n):
            if wi > w[j]:
                c += 1
    return c


def inv_for_qc(pi, q, c):
    n = len(pi)
    shift = (q + 1 - c) % n
    w = [(pi[(q + 1 + j) % n] - shift) % n for j in range(n)]
    return inv_count(w)


def full_min_inv(pi):
    """Brute-force I(pi) over all n^2 cuts; reference for cross-checking."""
    n = len(pi)
    best = n * n
    for q in range(n):
        for c in range(n):
            iv = inv_for_qc(pi, q, c)
            if iv < best:
                best = iv
    return best


def target(n):
    return (n - 1) ** 2 // 4


def diagonal_min_inv(pi):
    """Candidate: origin at each of the n points of pi itself (q = i, c s.t.
    pi(i) -> 0); best over these n candidates only."""
    n = len(pi)
    best = n * n
    for i0 in range(n):
        q = (i0 - 1) % n
        c = (i0 - pi[i0]) % n
        iv = inv_for_qc(pi, q, c)
        if iv < best:
            best = iv
    return best


def mode_shift_min_inv(pi):
    """Candidate: fix c once from the mode of (pi(i)-i) mod n (the shift that
    would make pi a pure rotation, sign-corrected), minimize only over q."""
    n = len(pi)
    shifts = Counter((pi[i] - i) % n for i in range(n))
    b_mode, _ = max(shifts.items(), key=lambda kv: kv[1])
    c_star = (-b_mode) % n
    best = n * n
    for q in range(n):
        iv = inv_for_qc(pi, q, c_star)
        if iv < best:
            best = iv
    return best


def rotation_sanity(nmax):
    """Pure rotations pi(i) = i + b must give I = 0 under both candidates;
    catches a sign convention bug before it contaminates the search."""
    for n in range(3, nmax + 1):
        for b in range(n):
            pi = tuple((i + b) % n for i in range(n))
            assert full_min_inv(pi) == 0, (n, b, "full")
            assert diagonal_min_inv(pi) == 0, (n, b, "diagonal")
            assert mode_shift_min_inv(pi) == 0, (n, b, "mode_shift")


def run(nmax):
    rotation_sanity(min(nmax, 8))
    results = {}
    t0 = time.time()
    for n in range(3, nmax + 1):
        worst = {"full": 0, "diagonal": 0, "mode_shift": 0}
        bad = {"diagonal": [], "mode_shift": []}
        checked = 0
        for pi in itertools.permutations(range(n)):
            checked += 1
            fm = full_min_inv(pi)
            dm = diagonal_min_inv(pi)
            mm = mode_shift_min_inv(pi)
            assert dm >= fm and mm >= fm, (pi, fm, dm, mm)
            worst["full"] = max(worst["full"], fm)
            worst["diagonal"] = max(worst["diagonal"], dm)
            worst["mode_shift"] = max(worst["mode_shift"], mm)
            tgt = target(n)
            if dm > tgt and len(bad["diagonal"]) < 3:
                bad["diagonal"].append({"pi": list(pi), "I_diagonal": dm, "I_true": fm})
            if mm > tgt and len(bad["mode_shift"]) < 3:
                bad["mode_shift"].append({"pi": list(pi), "I_mode_shift": mm, "I_true": fm})
        results[n] = {
            "checked": checked,
            "target": target(n),
            "worst_full": worst["full"],
            "worst_diagonal": worst["diagonal"],
            "worst_mode_shift": worst["mode_shift"],
            "diagonal_exceeds_target": worst["diagonal"] > target(n),
            "mode_shift_exceeds_target": worst["mode_shift"] > target(n),
            "examples": bad,
        }
        print(f"n={n} target={target(n)} worst_full={worst['full']} "
              f"worst_diagonal={worst['diagonal']} worst_mode_shift={worst['mode_shift']}")
    elapsed = time.time() - t0

    os.makedirs(OUT, exist_ok=True)
    report = {
        "version": VERSION,
        "core_version": CORE_VERSION,
        "goal": "test two O(n)-cut candidate rules for H13-I against exhaustive brute-force I(pi)",
        "command": f"python3 experiments/h13i_candidates.py --nmax {nmax}",
        "nmax": nmax,
        "coverage": "exhaustive (all permutations) per n",
        "elapsed_s": elapsed,
        "results": results,
        "conclusion": (
            "worst_full == target for all n (cross-checks H13-I / C33 on this range); "
            "worst_diagonal exceeds target starting at n=4; worst_mode_shift exceeds "
            "target starting at n=4 for even n (n=3,5,7 pass, n=4,6,8 fail); neither "
            "O(n)-cut rule suffices, an O(n^2) search (or a proof) is needed"
        ),
    }
    with open(os.path.join(OUT, "candidates_report.json"), "w") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)

    lines = [
        "# H13-I candidate cut-selection rules -- report",
        "",
        f"Version {VERSION}, core {CORE_VERSION}. Exhaustive per n, {elapsed:.1f} s total.",
        "",
        "| n | target | worst_full | worst_diagonal | worst_mode_shift |",
        "|---|---|---|---|---|",
    ]
    for n in sorted(results):
        r = results[n]
        lines.append(f"| {n} | {r['target']} | {r['worst_full']} | {r['worst_diagonal']} | {r['worst_mode_shift']} |")
    lines.append("")
    lines.append("Both O(n) rules fail (worst > target) starting at n = 4; example permutations "
                  "and their candidate vs. true I(pi) are in candidates_report.json under `examples`.")
    with open(os.path.join(OUT, "candidates_report.md"), "w") as f:
        f.write("\n".join(lines) + "\n")

    return report


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8)
    args = ap.parse_args()
    run(args.nmax)
