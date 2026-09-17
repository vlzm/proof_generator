"""Checker for C37 (toric classes, session 9): the averaging identity and its
corollaries, proved in docs/proofs/C37_toric_average.md.

Every assert below is tied to a numbered statement of that document.

  Part 1 (Lemma 1 + Lemma 7): for all pi with 4 <= n <= NMAX1 the table
    inv(a, s) computed from the definition (relabel the line, count inversions)
    agrees (i) with the pair criterion of Lemma 1 (count, for every pair, the
    cuts at which it is inverted, and compare the totals AB + (n-A)(n-B) /
    n(A+B) - 2AB) and (ii) with the incremental rule of Lemma 7 used by
    experiments/toric_min_inversions.c and experiments/toric_average.py.
  Part 2 (Lemma 2): the exact identity
    12 * sum_{a,s} inv(a,s) = 3 n^3 (n-1) + n^2 (n-1)(n-2) - 12 Delta(pi)
    for all pi with 4 <= n <= NMAX2 (integer arithmetic, no tolerance).
  Part 3 (Theorem 3 + Corollary 4): I(pi) <= floor((n-1)(2n-1)/6) for all pi
    with 4 <= n <= NMAX2, and the number of pi covered by Corollary 4
    (Delta(pi) >= n^2 (n^2 + 2 - 2*(n odd))/12 implies I(pi) <= floor((n-1)^2/4)):
    for those pi the conclusion is re-checked directly.
  Part 4 (Lemma 5 + Lemma 6): Delta(pi) = 0 exactly on the n reflections
    (4 <= n <= NMAX2); for 4 <= n <= NMAX4 and every h, inv(a,s) of pi_h equals
    C(K+1,2) + C(n-K-1,2) with K = (h-a-s) mod n at every cut, and
    I(pi_h) = floor((n-1)^2/4), attained at n cuts (n even) / 2n cuts (n odd).
  Part 5 (exhaustive scan, H13-I itself): the outputs of
    experiments/toric_min_inversions.c in data/runs/toric_inversions/ are
    re-read; max I = floor((n-1)^2/4) attained by exactly one toric class for
    4 <= n <= 13, and for 4 <= n <= NMAX2 the value of max I is recomputed here
    independently (definition, no incremental rule).  This part certifies the
    finite range of H13-I (C33), not the general statement.

Usage: python3 checks/check_C37.py [--nmax1 6] [--nmax2 8] [--nmax4 24]
Output: data/runs/check_C37/report.json, report.md.
"""

import argparse
import itertools
import json
import os
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
VERSION = "check_C37-1.0"
OUT = os.path.join(ROOT, "data", "runs", "check_C37")


def line(pi, a, s):
    n = len(pi)
    return [(pi[(a + j) % n] - s) % n for j in range(n)]


def inv_def(w):
    """Inversions of a linear sequence, straight from the definition."""
    return sum(1 for i in range(len(w)) for j in range(i + 1, len(w)) if w[i] > w[j])


def table_def(pi):
    n = len(pi)
    return [[inv_def(line(pi, a, s)) for s in range(n)] for a in range(n)]


def table_incr(pi):
    """Table by the incremental rule of Lemma 7."""
    n = len(pi)
    pinv = [0] * n
    for i, v in enumerate(pi):
        pinv[v] = i
    cur_a = inv_def(list(pi))
    tab = []
    for a in range(n):
        if a:
            cur_a += n - 1 - 2 * pi[a - 1]
        row, cur = [cur_a], cur_a
        for s in range(1, n):
            cur += n - 1 - 2 * ((pinv[s - 1] - a) % n)
            row.append(cur)
        tab.append(row)
    return tab


def defect(pi):
    n = len(pi)
    return sum(((i - j) % n + (pi[i] - pi[j]) % n - n) ** 2
               for i in range(n) for j in range(i + 1, n))


def part1(nmax, log):
    ok = True
    rows = []
    for n in range(4, nmax + 1):
        t0 = time.time()
        cnt = 0
        for pi in itertools.permutations(range(n)):
            td = table_def(pi)
            if td != table_incr(pi):
                ok = False
                log(f"FAIL Lemma 7: n={n} pi={pi}")
                break
            # Lemma 1: count per pair over all cuts
            tot_inv = sum(sum(r) for r in td)
            tot_pairs = 0
            for i in range(n):
                for j in range(i + 1, n):
                    A = (i - j) % n
                    B = (pi[i] - pi[j]) % n
                    tot_pairs += n * (A + B) - 2 * A * B
            if tot_inv != tot_pairs:
                ok = False
                log(f"FAIL Lemma 1: n={n} pi={pi} {tot_inv} != {tot_pairs}")
                break
            cnt += 1
        rows.append({"n": n, "perms": cnt, "seconds": round(time.time() - t0, 1)})
        log(f"part 1 (Lemmas 1, 7): n={n}, {cnt} permutations, definition = pair count = "
            f"incremental rule, {rows[-1]['seconds']} s")
    return ok, rows


def part2_3_4(nmax, log):
    ok = True
    rows = []
    for n in range(4, nmax + 1):
        t0 = time.time()
        target = (n - 1) ** 2 // 4
        weak = ((n - 1) * (2 * n - 1)) // 6
        thr = n * n * (n * n + 2 - 2 * (n % 2)) // 12   # Corollary 4
        covered = 0
        zero_defect = []
        cnt = 0
        for pi in itertools.permutations(range(n)):
            tab = table_incr(pi)
            tot = sum(sum(r) for r in tab)
            D = defect(pi)
            # Lemma 2
            if 12 * tot != 3 * n ** 3 * (n - 1) + n * n * (n - 1) * (n - 2) - 12 * D:
                ok = False
                log(f"FAIL Lemma 2: n={n} pi={pi}")
                break
            I = min(min(r) for r in tab)
            # Theorem 3
            if I > weak:
                ok = False
                log(f"FAIL Theorem 3: n={n} pi={pi} I={I} > {weak}")
                break
            # Corollary 4
            if D >= thr:
                covered += 1
                if I > target:
                    ok = False
                    log(f"FAIL Corollary 4: n={n} pi={pi} I={I} > {target}")
                    break
            if D == 0:
                zero_defect.append(pi)
            cnt += 1
        refl = [tuple((h - i) % n for i in range(n)) for h in range(n)]
        if sorted(zero_defect) != sorted(refl):
            ok = False
            log(f"FAIL Lemma 5: n={n} zero-defect set is not the set of reflections")
        rows.append({"n": n, "perms": cnt, "covered_by_corollary4": covered,
                     "threshold": thr, "weak_bound": weak, "target": target,
                     "zero_defect": len(zero_defect), "seconds": round(time.time() - t0, 1)})
        log(f"part 2-3 (Lemma 2, Theorem 3, Corollary 4): n={n}, {cnt} permutations, identity exact, "
            f"I <= {weak} (target {target}); Corollary 4 covers {covered}/{cnt}; "
            f"Lemma 5: {len(zero_defect)} permutations with Delta = 0 = the n reflections, "
            f"{rows[-1]['seconds']} s")
    return ok, rows


def part4_reflections(nmax, log):
    ok = True
    rows = []
    for n in range(4, nmax + 1):
        target = (n - 1) ** 2 // 4
        good_total = None
        for h in range(n):
            pi = tuple((h - i) % n for i in range(n))
            best, good = 10 ** 9, 0
            for a in range(n):
                for s in range(n):
                    K = (h - a - s) % n
                    got = inv_def(line(pi, a, s))
                    want = (K + 1) * K // 2 + (n - K - 1) * (n - K - 2) // 2
                    if got != want:
                        ok = False
                        log(f"FAIL Lemma 6 (formula): n={n} h={h} a={a} s={s} {got} != {want}")
                    if got < best:
                        best, good = got, 0
                    if got == best:
                        good += 1
            if best != target:
                ok = False
                log(f"FAIL Lemma 6: n={n} h={h} I={best} != {target}")
            exp_good = n if n % 2 == 0 else 2 * n
            if good != exp_good:
                ok = False
                log(f"FAIL Lemma 6 (count): n={n} h={h} {good} != {exp_good}")
            good_total = good
        rows.append({"n": n, "I_reflection": target, "optimal_cuts": good_total})
        log(f"part 4 (Lemma 6): n={n}, all h: inv = C(K+1,2)+C(n-K-1,2) at all {n * n} cuts, "
            f"I(pi_h) = {target} = floor((n-1)^2/4) at {good_total} cuts")
    return ok, rows


def part5(nmax_indep, log):
    ok = True
    rows = []
    for n in range(4, 14):
        path = os.path.join(ROOT, "data", "runs", "toric_inversions", f"toric_n{n}.json")
        if not os.path.exists(path):
            continue
        r = json.load(open(path))
        target = (n - 1) ** 2 // 4
        good = (r["max_I"] == target and r["max_I_count"] == 1 and r["target"] == target)
        indep = None
        if n <= nmax_indep:
            indep = max(min(min(row) for row in table_def(pi))
                        for pi in itertools.permutations(range(n)) if pi[0] == 0)
            good = good and indep == r["max_I"]
        ok = ok and good
        rows.append({"n": n, "reps": r["reps"], "max_I": r["max_I"], "target": target,
                     "classes_at_max": r["max_I_count"], "min_good_cuts": r["min_good_cuts"],
                     "independent_max_I": indep, "seconds": r["seconds"]})
        log(f"part 5 (H13-I finite range): n={n}, {r['reps']} class representatives, max I = {r['max_I']} "
            f"= floor((n-1)^2/4) = {target}, attained by {r['max_I_count']} toric class"
            f"{'' if indep is None else f', independent recomputation {indep}'}"
            f"; min number of cuts with inv <= target = {r['min_good_cuts']} -> {'ok' if good else 'FAIL'}")
    return ok, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax1", type=int, default=6)
    ap.add_argument("--nmax2", type=int, default=8)
    ap.add_argument("--nmax4", type=int, default=24)
    ap.add_argument("--nmax5", type=int, default=7)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    lines = []

    def log(msg):
        print(msg, flush=True)
        lines.append(msg)

    log(f"== {VERSION} args={vars(args)}")
    ok1, r1 = part1(args.nmax1, log)
    ok2, r2 = part2_3_4(args.nmax2, log)
    ok4, r4 = part4_reflections(args.nmax4, log)
    ok5, r5 = part5(args.nmax5, log)
    verdict = "PASS" if (ok1 and ok2 and ok4 and ok5) else "FAIL"
    log(f"verdict: {verdict}")
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump({"version": VERSION, "args": vars(args), "part1": r1, "part2_3": r2,
                   "part4": r4, "part5": r5, "verdict": verdict}, f, indent=1)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# check_C37 — среднее по разрезам торического класса (леммы 1–7, теорема 3)\n\n")
        f.write("Команда: `python3 checks/check_C37.py --nmax1 %d --nmax2 %d --nmax4 %d --nmax5 %d`. "
                "Версия: %s.\nДоказательство: `docs/proofs/C37_toric_average.md`.\n\n```text\n"
                % (args.nmax1, args.nmax2, args.nmax4, args.nmax5, VERSION))
        f.write("\n".join(lines) + "\n```\n")
    sys.exit(0 if verdict == "PASS" else 1)


if __name__ == "__main__":
    main()
