"""Independent checker for C37 and C38 (docs/proofs/C37_toric_identity.md).

Everything here is recomputed from the definitions: inversions of a line are
counted by the naive double loop over pairs, D2 by the definition, and the cut
tables by re-linearising the permutation for each of the n^2 cuts.  The step
recurrences of lemmas 1-2 are *checked*, never used to build the tables.

  part 1 (lemma 1, lemma 2): for every pi and every cut,
        inv(a+1,b) - inv(a,b) = n-1-2 w_0,  D2(a+1,b) - D2(a,b) = n(n-1-2 w_0),
        inv(a,b+1) - inv(a,b) = n-1-2 p,    D2(a,b+1) - D2(a,b) = n(n-1-2 p),
      where w_0 is the first letter of the line and p the position of letter 0.
  part 2 (theorem A): inv(a,b) - D2(a,b)/n is one and the same rational for all
      n^2 cuts, and equals (n^2-1)/12 - W/n^2 with W the sawtooth correlation.
  part 3 (lemma 3): for every unordered pair {i,j} the number of cuts inverting
      it equals n(d+e) - 2de, counted by brute force over all n^2 cuts.
  part 4 (theorem B): max over S_n of A(pi) equals (n-1)(2n-1)/6 and the set of
      maximisers is exactly the n reflections; scalar identity
      sum_{k=1}^{n-1}(k-n/2)^2 = n(n-1)(n-2)/12 in exact arithmetic, 2 <= n <= 200.
  part 5 (corollary B1): I(pi) <= floor((n-1)(2n-1)/6) for every pi (and the
      conjectured H13-I value floor((n-1)^2/4) is printed alongside).

Usage: python3 checks/check_C37.py [--full 7] [--exhaustive 8] [--sample 300]
  --full        : parts 1-3 exhaustively over S_n for 4 <= n <= full
  --exhaustive  : parts 4-5 exhaustively over S_n for 4 <= n <= exhaustive
  --sample      : parts 1-2 on this many random pi for n = full+1, full+2
Output: data/runs/check_C37/report.json, report.md.
"""

import argparse
import itertools
import json
import os
import random
import time
from fractions import Fraction

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
VERSION = "check_C37-1.0"
OUT = os.path.join(ROOT, "data", "runs", "check_C37")


def line(pi, a, b):
    n = len(pi)
    return [(pi[(j + a) % n] - b) % n for j in range(n)]


def inv(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def d2(w):
    return sum((j - w[j]) ** 2 for j in range(len(w)))


def sawtooth_const(pi):
    n = len(pi)
    f = [Fraction(2 * t - (n - 1), 2) for t in range(n)]
    W = sum(f[(j - i) % n] * f[(pi[j] - pi[i]) % n] for i in range(n) for j in range(n))
    return Fraction(n * n - 1, 12) - Fraction(W, n * n)


def is_reflection(pi):
    n = len(pi)
    s = pi[0] % n
    return all((pi[i] + i) % n == s for i in range(n))


def part12(pi, log):
    """Lemmas 1, 2 and theorem A for one permutation.  Returns the constant."""
    n = len(pi)
    iv = {}
    dq = {}
    for a in range(n):
        for b in range(n):
            w = line(pi, a, b)
            iv[(a, b)] = inv(w)
            dq[(a, b)] = d2(w)
    consts = set()
    for a in range(n):
        for b in range(n):
            w = line(pi, a, b)
            w0 = w[0]
            p = w.index(0)
            assert iv[((a + 1) % n, b)] - iv[(a, b)] == n - 1 - 2 * w0, ("L1 inv", pi, a, b)
            assert dq[((a + 1) % n, b)] - dq[(a, b)] == n * (n - 1 - 2 * w0), ("L1 D2", pi, a, b)
            assert iv[(a, (b + 1) % n)] - iv[(a, b)] == n - 1 - 2 * p, ("L2 inv", pi, a, b)
            assert dq[(a, (b + 1) % n)] - dq[(a, b)] == n * (n - 1 - 2 * p), ("L2 D2", pi, a, b)
            consts.add(Fraction(iv[(a, b)]) - Fraction(dq[(a, b)], n))
    assert len(consts) == 1, ("theorem A", pi, sorted(consts)[:4])
    c = consts.pop()
    assert c == sawtooth_const(pi), ("const formula", pi, c, sawtooth_const(pi))
    return c, iv


def part3(pi, iv):
    """Lemma 3: number of cuts inverting a given pair."""
    n = len(pi)
    for i in range(n):
        for j in range(i + 1, n):
            d = (j - i) % n
            e = (pi[j] - pi[i]) % n
            cnt = 0
            for a in range(n):
                for b in range(n):
                    xi, xj = (i - a) % n, (j - a) % n
                    yi, yj = (pi[i] - b) % n, (pi[j] - b) % n
                    if (xi < xj) != (yi < yj):
                        cnt += 1
            assert cnt == n * (d + e) - 2 * d * e, ("L3", pi, i, j, cnt)
    # consistency: the pair counts add up to sum of inv over all cuts
    tot = sum(n * (((j - i) % n) + ((pi[j] - pi[i]) % n))
              - 2 * ((j - i) % n) * ((pi[j] - pi[i]) % n)
              for i in range(n) for j in range(i + 1, n))
    assert tot == sum(iv.values()), ("L4", pi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", type=int, default=7)
    ap.add_argument("--exhaustive", type=int, default=8)
    ap.add_argument("--sample", type=int, default=300)
    ap.add_argument("--seed", type=int, default=20260913)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    rng = random.Random(args.seed)
    log = []
    t00 = time.time()
    log.append(f"== {VERSION} args={vars(args)}")

    # parts 1-3, exhaustive
    for n in range(4, args.full + 1):
        t0 = time.time()
        cnt = 0
        for pi in itertools.permutations(range(n)):
            c, iv = part12(pi, log)
            part3(pi, iv)
            cnt += 1
        log.append(f"parts 1-3 n={n}: {cnt} pi x {n*n} разрезов, леммы 1-3, теорема A "
                   f"и формула константы — ok, {time.time()-t0:.1f} s")

    # parts 1-2, sampled for two larger n
    for n in range(args.full + 1, args.full + 3):
        t0 = time.time()
        for _ in range(args.sample):
            pi = list(range(n))
            rng.shuffle(pi)
            part12(tuple(pi), log)
        log.append(f"parts 1-2 n={n}: {args.sample} случайных pi (seed {args.seed}) — ok, "
                   f"{time.time()-t0:.1f} s")

    # part 4: scalar identity
    for n in range(2, 201):
        lhs = sum(Fraction(2 * k - n, 2) ** 2 for k in range(1, n))
        assert lhs == Fraction(n * (n - 1) * (n - 2), 12), ("scalar", n)
    log.append("part 4: sum_{k=1}^{n-1}(k-n/2)^2 = n(n-1)(n-2)/12 точно, 2 <= n <= 200 — ok")

    # parts 4-5: exhaustive maximum of the mean and corollary B1
    rows = []
    for n in range(4, args.exhaustive + 1):
        t0 = time.time()
        best = -1
        argmax = []
        maxI = -1
        bound = (n - 1) * (2 * n - 1) // 6
        for pi in itertools.permutations(range(n)):
            tot = 0
            I = None
            for a in range(n):
                for b in range(n):
                    v = inv(line(pi, a, b))
                    tot += v
                    if I is None or v < I:
                        I = v
            assert I <= bound, ("B1", pi, I, bound)
            maxI = max(maxI, I)
            if tot > best:
                best, argmax = tot, [pi]
            elif tot == best:
                argmax.append(pi)
        exp = n * n * (n - 1) * (2 * n - 1) // 6
        assert best == exp, ("theorem B", n, best, exp)
        assert len(argmax) == n and all(is_reflection(p) for p in argmax), ("equality", n)
        rows.append({"n": n, "max_n2A": best, "bound_n2": exp, "argmax": len(argmax),
                     "maxI": maxI, "B1": bound, "M": ((n - 1) ** 2) // 4})
        log.append(f"parts 4-5 n={n}: max n^2 A = {best} = n^2 (n-1)(2n-1)/6, "
                   f"argmax = {len(argmax)} отражений; max I = {maxI} <= "
                   f"floor((n-1)(2n-1)/6) = {bound} (H13-I: {((n-1)**2)//4}), "
                   f"{time.time()-t0:.1f} s")

    log.append("verdict: PASS")
    text = "\n".join(log)
    print(text)
    with open(os.path.join(OUT, "report.json"), "w") as fh:
        json.dump({"version": VERSION, "args": vars(args), "rows": rows,
                   "log": log, "seconds": round(time.time() - t00, 1)}, fh, indent=1)
    with open(os.path.join(OUT, "report.md"), "w") as fh:
        fh.write("# check_C37 — тождество Кендалл–Спирмен и максимум среднего\n\n"
                 f"Команда: `python3 checks/check_C37.py --full {args.full} "
                 f"--exhaustive {args.exhaustive} --sample {args.sample}`. "
                 f"Версии: {VERSION}. Всё пересчитано от определений "
                 "(инверсии — перебором пар).\n\n```text\n" + text +
                 f"\n```\n\nВремя: {round(time.time()-t00, 1)} с.\n")


if __name__ == "__main__":
    main()
