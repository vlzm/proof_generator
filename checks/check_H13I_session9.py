"""Checker for the session-9 work on H13-I (docs/notes/h13_line_model.md §7).

H13-I is NOT proved this session; this checker verifies only what was
established: the new fast exhaustive implementation, the extended range,
and the two new small lemmas (pair-weight identity, inverse symmetry).
It does not certify H13-I itself (rule 5: exhaustive checking is not a proof
of the general case).

  Part 1 (independent re-implementation, repo rule 3): compile both
    experiments/line_profile.c (slow O(n^4)/permutation reference, session 8)
    and experiments/line_profile_fast.c (O(n^2)/permutation, session 9); run
    both for 4 <= n <= AMAX and compare the full histogram of I(pi) values
    (count per I) exactly.
  Part 2 (Lemma D, pair-weight identity): for random pi at several n, check
    sum_{(q,c)} inv(w_{q,c}) == sum_{i<j} f(dx(i,j), dy(i,j)) with
    f(dx,dy) = n*(dx+dy) - 2*dx*dy, dx=(j-i) mod n, dy=(pi(j)-pi(i)) mod n,
    against a brute-force sum over all n^2 cuts.
  Part 3 (Lemma E, I(pi) = I(pi^{-1})): exhaustive for 3 <= n <= AMAX_INV.
  Part 4 (ruled-out one-parameter families, for the record): exhaustively
    confirm the known violations of
      G(pi) = min_q inv(pi read from q), no value shift,
      J(pi) = min_c inv(pi - c mod n), no position rotation,
    against the fixed reference numbers recorded in the note (regression
    check on the counterexamples, not a claim that G or J are bounded).
  Part 5 (extended exhaustive range of C33): read
    data/runs/line_profile_fast/hist_n{n}.json for n up to whatever was
    computed (4 <= n <= 12, plus n = 13 if present) and check
    max_I == floor((n-1)^2/4).

Usage: python3 checks/check_H13I_session9.py [--amax 9]
Output: data/runs/check_H13I_session9/report.json, report.md.
"""

import argparse
import itertools
import json
import math
import os
import random
import subprocess
import sys
import tempfile
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "data", "runs", "check_H13I_session9")
VERSION = "check_H13I_session9-1.0"


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def I_of(pi):
    n = len(pi)
    best = None
    for q in range(n):
        for c in range(n):
            shift = q + 1 - c
            w = [(pi[(q + 1 + j) % n] - shift) % n for j in range(n)]
            iv = inversions(w)
            if best is None or iv < best:
                best = iv
    return best


def part1(amax, log, tmpdir):
    """Cross-check the fast implementation against the slow session-8 reference."""
    slow_bin = os.path.join(tmpdir, "line_profile_slow")
    fast_bin = os.path.join(tmpdir, "line_profile_fast")
    subprocess.run(["gcc", "-O2", "-o", slow_bin, os.path.join(ROOT, "experiments", "line_profile.c")], check=True)
    subprocess.run(["gcc", "-O2", "-o", fast_bin, os.path.join(ROOT, "experiments", "line_profile_fast.c")], check=True)
    ok = True
    rows = []
    for n in range(4, amax + 1):
        total = math.factorial(n)
        dist_path = os.path.join(tmpdir, f"zero_dist_n{n}.bin")
        with open(dist_path, "wb") as f:
            f.write(bytes(total))
        slow = json.loads(subprocess.run([slow_bin, str(n), dist_path], capture_output=True, text=True, check=True).stdout)
        fast = json.loads(subprocess.run([fast_bin, str(n)], capture_output=True, text=True, check=True).stdout)
        os.remove(dist_path)
        slow_hist = {row["I"]: row["count"] for row in slow["by_I"]}
        fast_hist = dict(fast["hist"])
        same = slow_hist == fast_hist
        ok = ok and same
        rows.append({"n": n, "match": same, "max_I_slow": max(slow_hist), "max_I_fast": fast["max_I"]})
        log(f"part 1 n={n}: histograms {'match' if same else 'MISMATCH'} "
            f"(max I slow={max(slow_hist)} fast={fast['max_I']})")
    return ok, rows


def part2(log, seed=0):
    """Lemma D: pair-weight identity for the sum over all n^2 cuts."""
    random.seed(seed)
    ok = True
    rows = []
    for n in range(3, 9):
        for _ in range(5):
            pi = list(range(n))
            random.shuffle(pi)
            brute = 0
            for q in range(n):
                for c in range(n):
                    shift = q + 1 - c
                    w = [(pi[(q + 1 + j) % n] - shift) % n for j in range(n)]
                    brute += inversions(w)
            predicted = 0
            for i in range(n):
                for j in range(i + 1, n):
                    dx = (j - i) % n
                    dy = (pi[j] - pi[i]) % n
                    predicted += n * (dx + dy) - 2 * dx * dy
            good = brute == predicted
            ok = ok and good
            if not good:
                rows.append({"n": n, "pi": pi, "brute": brute, "predicted": predicted, "ok": good})
    log(f"part 2 (Lemma D, pair-weight identity): {'PASS' if ok else 'FAIL'} "
        f"(5 random pi for each 3 <= n <= 8, brute sum over all n^2 cuts vs closed form)")
    return ok, rows


def part3(amax_inv, log):
    """Lemma E: I(pi) = I(pi^{-1})."""
    ok = True
    rows = []
    for n in range(3, amax_inv + 1):
        bad = 0
        checked = 0
        for perm in itertools.permutations(range(n)):
            pi = list(perm)
            inv_pi = [0] * n
            for i, v in enumerate(pi):
                inv_pi[v] = i
            if I_of(pi) != I_of(inv_pi):
                bad += 1
            checked += 1
        good = bad == 0
        ok = ok and good
        rows.append({"n": n, "checked": checked, "mismatches": bad})
        log(f"part 3 (Lemma E, I(pi)=I(pi^-1)) n={n}: {checked} pi, {bad} mismatches -> {'PASS' if good else 'FAIL'}")
    return ok, rows


REFERENCE_G = {4: 2, 5: 4, 6: 6, 7: 10, 8: 13, 9: 17}
REFERENCE_J = {4: 2, 5: 4, 6: 6, 7: 10, 8: 13, 9: 17}


def part4(amax, log):
    """Ruled-out one-parameter families G (position rotation only) and J
    (value shift only): regression-check the known worst-case values against
    the ones recorded in the note (both first exceed the bound at n = 7)."""
    ok = True
    rows = []
    for n in range(4, amax + 1):
        maxG = -1
        maxJ = -1
        for perm in itertools.permutations(range(n)):
            pi = list(perm)
            g = min(inversions(pi[q:] + pi[:q]) for q in range(n))
            j = min(inversions([(x - c) % n for x in pi]) for c in range(n))
            maxG = max(maxG, g)
            maxJ = max(maxJ, j)
        bound = (n - 1) ** 2 // 4
        good = (maxG == REFERENCE_G.get(n, maxG)) and (maxJ == REFERENCE_J.get(n, maxJ))
        ok = ok and good
        rows.append({"n": n, "max_G": maxG, "max_J": maxJ, "bound": bound,
                     "G_exceeds_bound": maxG > bound, "J_exceeds_bound": maxJ > bound})
        log(f"part 4 n={n}: max G(pi) (position-rotation only) = {maxG}, "
            f"max J(pi) (value-shift only) = {maxJ}, bound = {bound} "
            f"({'exceeds' if maxG > bound else 'within'} / {'exceeds' if maxJ > bound else 'within'}) "
            f"-> {'as recorded' if good else 'CHANGED from recorded values'}")
    return ok, rows


def part5(log):
    """Read the extended exhaustive range (session 9, n up to 12 or 13) and
    check max_I == floor((n-1)^2/4)."""
    ok = True
    rows = []
    d = os.path.join(ROOT, "data", "runs", "line_profile_fast")
    for n in range(4, 20):
        path = os.path.join(d, f"hist_n{n}.json")
        if not os.path.exists(path):
            continue
        r = json.load(open(path))
        bound = (n - 1) ** 2 // 4
        good = r["max_I"] == bound
        ok = ok and good
        rows.append({"n": n, "count": r["count"], "max_I": r["max_I"], "bound": bound, "ok": good})
        log(f"part 5 n={n}: max I(pi) = {r['max_I']} (floor((n-1)^2/4) = {bound}), "
            f"{r['count']} pi -> {'ok' if good else 'FAIL'}")
    return ok, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amax", type=int, default=9, help="max n for part 1 (slow/fast cross-check)")
    ap.add_argument("--amax_inv", type=int, default=7, help="max n for part 3 (I(pi)=I(pi^-1))")
    ap.add_argument("--amax_g", type=int, default=9, help="max n for part 4 (G, J families)")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    lines = []

    def log(msg):
        print(msg, flush=True)
        lines.append(msg)

    log(f"== {VERSION} args={vars(args)}")
    t0 = time.time()
    with tempfile.TemporaryDirectory() as tmpdir:
        ok1, r1 = part1(args.amax, log, tmpdir)
    ok2, r2 = part2(log)
    ok3, r3 = part3(args.amax_inv, log)
    ok4, r4 = part4(args.amax_g, log)
    ok5, r5 = part5(log)
    verdict = "PASS" if (ok1 and ok2 and ok3 and ok4 and ok5) else "FAIL"
    log(f"total time: {time.time() - t0:.1f} s")
    log(f"verdict: {verdict}  (note: this PASS certifies the session-9 lemmas and computational "
        f"data, NOT H13-I itself, which remains unproven -- see h13_line_model.md par 7)")
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump({"version": VERSION, "args": vars(args), "part1_fast_vs_slow": r1,
                   "part2_pair_weight_identity": r2, "part3_inverse_symmetry": r3,
                   "part4_ruled_out_families": r4, "part5_extended_range": r5,
                   "verdict": verdict}, f, indent=1)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# check_H13I_session9 — новые леммы и расширенный диапазон (H13-I не доказана)\n\n")
        f.write("Команда: `python3 checks/check_H13I_session9.py --amax %d`. Версия: %s.\n\n```text\n" %
                (args.amax, VERSION))
        f.write("\n".join(lines) + "\n```\n")
    sys.exit(0 if verdict == "PASS" else 1)


if __name__ == "__main__":
    main()
