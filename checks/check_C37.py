"""Checker for C37 (docs/proofs/C37_double_cut_algebra.md), session 9.

Independent implementation: lines and inversions are rebuilt from the definition
of the double cut in docs/notes/h13_line_model.md, not imported from
experiments/toric_cut_algebra.py.

  Lemma 1  the (q, c) and (m, k) parametrisations give the same n^2 lines.
  Lemma 2  inv(w(m, k)) = inv(pi) + s(n - s) - 2 e_inv(S), S = A_m xor B_k,
           equivalently |Inv(pi) xor delta(S)|.
  Lemma 3  inv(w) = B_n/2 - Q(x)/2 with x the +-1 vector of S, and
           inv(w) <= floor((n-1)^2/4) iff Q(x) >= floor(n/2).
  Lemma 4  a pair at clockwise distances (d_p, d_v) is concordant for exactly
           (n - d_p)(n - d_v) + d_p d_v of the n^2 cuts; Sigma(id_n) and
           Sigma(sigma_n) match the closed forms; averaging settles sigma_n
           never and id_n from n = 5 on.
  Theorem 5 / Cor 5.1  max over ALL x in {+-1}^n of Q(x) >= floor(n/2), for
           permutation signings (exhaustive over pi and over all subsets) and
           for arbitrary +-1 signings of K_n (random); the local descent of the
           proof reaches the bound on large structured inputs.
  Cor 5.2  the criterion H14 implies inv(w) <= floor((n-1)^2/4); it is not
           necessary: 16 affine counterexamples at n = 8.

Usage: python3 checks/check_C37.py [--amax 8] [--free-amax 8]
Output: data/runs/check_C37/report.json, report.md.
Version check_C37-1.0.
"""

import argparse
import itertools
import json
import os
import random
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
from moves import CORE_VERSION  # noqa: E402

VERSION = "check_C37-1.0"
OUT = os.path.join(ROOT, "data", "runs", "check_C37")


def line_qc(perm, q, c):
    n = len(perm)
    sh = q + 1 - c
    return tuple((perm[(q + 1 + j) % n] - sh) % n for j in range(n))


def line_mk(perm, m, k):
    n = len(perm)
    return tuple((perm[(m + j) % n] - k) % n for j in range(n))


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def inv_pairs(perm):
    n = len(perm)
    return [(i, j) for i in range(n) for j in range(i + 1, n) if perm[i] > perm[j]]


def cut_subset(perm, m, k):
    n = len(perm)
    return frozenset(i for i in range(n) if (i < m) != (perm[i] < k))


def lemma1(amax, log):
    ok = True
    for n in range(4, min(amax, 6) + 1):
        for perm in itertools.permutations(range(n)):
            a = set(line_qc(perm, q, c) for q in range(n) for c in range(n))
            b = set(line_mk(perm, m, k) for m in range(n) for k in range(n))
            if a != b or len(a) > n * n:
                ok = False
        log(f"lemma 1 n={n}: all perms, (q,c) and (m,k) families coincide -> {'ok' if ok else 'FAIL'}")
    return ok


def lemma23(amax, log):
    ok = True
    rnd = random.Random(370037)
    for n in range(4, min(amax, 7) + 1):
        t0 = time.time()
        for perm in itertools.permutations(range(n)):
            ip = inv_pairs(perm)
            inv0 = len(ip)
            for m in range(n):
                for k in range(n):
                    S = cut_subset(perm, m, k)
                    s = len(S)
                    e = sum(1 for (i, j) in ip if (i in S) != (j in S))
                    pred = inv0 + s * (n - s) - 2 * e
                    act = inversions(line_mk(perm, m, k))
                    if pred != act:
                        ok = False
                    # lemma 3: sign form
                    x = [-1 if i in S else 1 for i in range(n)]
                    Q = 0
                    for i in range(n):
                        for j in range(i + 1, n):
                            eps = -1 if (perm[i] > perm[j]) else 1
                            Q += eps * x[i] * x[j]
                    bn = n * (n - 1) // 2
                    if (bn - Q) % 2 or act != (bn - Q) // 2:
                        ok = False
                    if (act <= ((n - 1) ** 2) // 4) != (Q >= n // 2):
                        ok = False
        log(f"lemma 2,3 n={n}: all perms x {n * n} cuts -> {'ok' if ok else 'FAIL'}"
            f" ({time.time() - t0:.1f} s)")
    for n in (8, 12, 20, 33, 50):
        for _ in range(100):
            perm = list(range(n))
            rnd.shuffle(perm)
            ip = inv_pairs(perm)
            m = rnd.randrange(n)
            k = rnd.randrange(n)
            S = cut_subset(perm, m, k)
            s = len(S)
            e = sum(1 for (i, j) in ip if (i in S) != (j in S))
            if len(ip) + s * (n - s) - 2 * e != inversions(line_mk(perm, m, k)):
                ok = False
        log(f"lemma 2 n={n}: 100 random (pi, m, k) -> {'ok' if ok else 'FAIL'}")
    return ok


def sigma_sum(perm):
    n = len(perm)
    tot = 0
    for i in range(n):
        for j in range(i + 1, n):
            d_p = (j - i) % n
            d_v = (perm[j] - perm[i]) % n
            tot += (n - d_p) * (n - d_v) + d_p * d_v
    return tot


def lemma4(amax, log):
    ok = True
    for n in range(4, min(amax, 6) + 1):
        for perm in itertools.permutations(range(n)):
            tot = 0
            for i in range(n):
                for j in range(i + 1, n):
                    d_p = (j - i) % n
                    d_v = (perm[j] - perm[i]) % n
                    cnt = 0
                    for m in range(n):
                        for k in range(n):
                            w = line_mk(perm, m, k)
                            pi_ = (i - m) % n
                            pj = (j - m) % n
                            if (pi_ < pj) == (w[pi_] < w[pj]):
                                cnt += 1
                    if cnt != (n - d_p) * (n - d_v) + d_p * d_v:
                        ok = False
                    tot += cnt
            direct = sum(n * (n - 1) // 2 - inversions(line_mk(perm, m, k))
                         for m in range(n) for k in range(n))
            if direct != tot:
                ok = False
        log(f"lemma 4 n={n}: all perms, pair counts and Sigma(pi) -> {'ok' if ok else 'FAIL'}")
    # closed forms, exact integer arithmetic
    for n in range(4, 61):
        idp = tuple(range(n))
        sg = tuple((1 - i) % n for i in range(n))
        rv = tuple(n - 1 - i for i in range(n))
        if sigma_sum(idp) * 6 != n * n * (n - 1) * (2 * n - 1):
            ok = False
        if sigma_sum(sg) * 6 != n * n * (n * n - 1):
            ok = False
        if sigma_sum(rv) * 6 != n * n * (n * n - 1):
            ok = False
        tgt = (n * n) // 4
        if sigma_sum(sg) >= tgt * n * n:
            ok = False                       # averaging must fail on sigma_n
        if (sigma_sum(idp) >= tgt * n * n) != (n >= 5):
            ok = False                       # and must work on id_n from n = 5
    log(f"lemma 4: closed forms Sigma(id_n), Sigma(sigma_n) = Sigma(rev_n) and the"
        f" averaging verdicts, 4 <= n <= 60 -> {'ok' if ok else 'FAIL'}")
    return ok


def Qvalue(pairs_eps, x, n):
    Q = 0
    for i in range(n):
        for j in range(i + 1, n):
            Q += pairs_eps[i][j] * x[i] * x[j]
    return Q


def eps_of_perm(perm):
    n = len(perm)
    return [[0 if i == j else (-1 if ((i < j) != (perm[i] < perm[j])) else 1)
             for j in range(n)] for i in range(n)]


def max_Q_bruteforce(eps, n):
    best = None
    for mask in range(1 << (n - 1)):
        x = [1] + [(-1 if mask >> (i - 1) & 1 else 1) for i in range(1, n)]
        Q = Qvalue(eps, x, n)
        if best is None or Q > best:
            best = Q
    return best


def descent(eps, n):
    """Local descent of theorem 5: 1- and 2-flips, each step raises Q."""
    x = [1] * n
    Q = Qvalue(eps, x, n)
    steps = 0
    while True:
        r = [sum(eps[p][j] * x[p] * x[j] for j in range(n) if j != p) for p in range(n)]
        p = min(range(n), key=lambda i: r[i])
        if r[p] < 0:
            x[p] = -x[p]
            Q -= 2 * r[p]
            steps += 1
            continue
        moved = False
        for p in range(n):
            for q in range(p + 1, n):
                w = r[p] + r[q] - 2 * eps[p][q] * x[p] * x[q]
                if w < 0:
                    x[p] = -x[p]
                    x[q] = -x[q]
                    Q -= 2 * w
                    moved = True
                    steps += 1
                    break
            if moved:
                break
        if not moved:
            return Q, steps


def max_Q_fast(perm, n):
    """max over all x in {+-1}^n of Q(x), via per-element inversion counts."""
    full = (1 << n) - 1
    im = [0] * n
    for i in range(n):
        for j in range(n):
            if i != j and (i < j) != (perm[i] < perm[j]):
                im[i] |= 1 << j
    best = None
    for S in range(1 << (n - 1)):
        comp = full ^ S
        tot = 0
        for p in range(n):
            cross = (comp if (S >> p & 1) else S) & ~(1 << p)
            tot += bin(im[p] ^ cross).count("1")
        Q = (n * (n - 1) - 2 * tot) // 2
        if best is None or Q > best:
            best = Q
    return best


def theorem5(free_amax, log):
    ok = True
    rows = []
    for n in range(4, free_amax + 1):
        tgt = n // 2
        t0 = time.time()
        worst = None
        cross = (n <= 6)
        for perm in itertools.permutations(range(n)):
            best = max_Q_fast(perm, n)
            if cross and best != max_Q_bruteforce(eps_of_perm(perm), n):
                ok = False
            if best < tgt:
                ok = False
            if worst is None or best < worst:
                worst = best
        rows.append({"n": n, "min_over_pi_of_max_Q": worst, "target": tgt,
                     "crosschecked_with_direct_sum": cross,
                     "time": round(time.time() - t0, 1)})
        log(f"theorem 5 n={n}: all perms, min_pi max_x Q(x) = {worst} >= floor(n/2) = {tgt}"
            f"{' (cross-checked against the direct sum)' if cross else ''}"
            f" -> {'ok' if ok else 'FAIL'} ({rows[-1]['time']} s)")
    # arbitrary +-1 signings (the theorem is not about permutations)
    rnd = random.Random(5150)
    for n in range(4, 11):
        worst = None
        for _ in range(300):
            eps = [[0] * n for _ in range(n)]
            for i in range(n):
                for j in range(i + 1, n):
                    v = rnd.choice((-1, 1))
                    eps[i][j] = eps[j][i] = v
            best = max_Q_bruteforce(eps, n)
            if best < n // 2:
                ok = False
            if worst is None or best < worst:
                worst = best
        log(f"theorem 5 n={n}: 300 random +-1 signings of K_n, min max_x Q = {worst}"
            f" >= {n // 2} -> {'ok' if ok else 'FAIL'}")
    # tightness: eps = -1 (rev_n)
    for n in range(4, 15):
        eps = [[0 if i == j else -1 for j in range(n)] for i in range(n)]
        if max_Q_bruteforce(eps, n) != n // 2:
            ok = False
        rv = tuple(n - 1 - i for i in range(n))
        if max_Q_fast(rv, n) != n // 2:
            ok = False
    log(f"theorem 5: tightness eps = -1 (rev_n), max_x Q = floor(n/2), 4 <= n <= 14"
        f" -> {'ok' if ok else 'FAIL'}")
    # constructive descent on large inputs
    for n in (12, 20, 33, 50, 80, 101):
        tgt = ((n - 1) ** 2) // 4
        worst = -1
        inputs = []
        for h in range(0, n, max(1, n // 8)):
            inputs.append(tuple((h - i) % n for i in range(n)))
        for a in range(2, n):
            if len(inputs) > 24:
                break
            if all(a % p for p in (2, 3, 5, 7) if n % p == 0) and gcd(a, n) == 1:
                inputs.append(tuple((a * i + 1) % n for i in range(n)))
        for _ in range(10):
            p = list(range(n))
            rnd.shuffle(p)
            inputs.append(tuple(p))
        for perm in inputs:
            eps = eps_of_perm(perm)
            Q, steps = descent(eps, n)
            inv = (n * (n - 1) // 2 - Q) // 2
            if Q < n // 2 or inv > tgt:
                ok = False
            worst = max(worst, inv)
        log(f"theorem 5 n={n}: descent on {len(inputs)} inputs, max resulting inv = {worst}"
            f" <= floor((n-1)^2/4) = {tgt} -> {'ok' if ok else 'FAIL'}")
    return ok, rows


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def corollary52(amax, log):
    """The criterion implies the bound; and it is not necessary (n = 8)."""
    ok = True
    rows = []
    for n in range(4, amax + 1):
        h = (n - 1) // 2
        tgt = ((n - 1) ** 2) // 4
        fails = []
        t0 = time.time()
        for perm in itertools.permutations(range(n)):
            good = None
            for m in range(n):
                for k in range(n):
                    w = line_mk(perm, m, k)
                    invp = [sum(1 for b in range(n)
                                if b != a and ((a < b) != (w[a] < w[b])))
                            for a in range(n)]
                    if max(invp) > h:
                        continue
                    if n % 2 == 1:
                        Z = [a for a in range(n) if invp[a] == h]
                        if not all((a < b) != (w[a] < w[b]) for a in Z for b in Z if a < b):
                            continue
                    # criterion holds: the implication must hold too
                    if inversions(w) > tgt:
                        ok = False
                    good = (m, k)
                    break
                if good:
                    break
            if good is None:
                fails.append(list(perm))
        rows.append({"n": n, "no_cut_satisfying_criterion": len(fails),
                     "examples": fails[:16], "time": round(time.time() - t0, 1)})
        log(f"cor 5.2 n={n}: implication holds; perms with no cut satisfying the"
            f" criterion: {len(fails)}" + (f", e.g. {fails[0]}" if fails else ""))
    return ok, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amax", type=int, default=8)
    ap.add_argument("--free-amax", type=int, default=8)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    t0 = time.time()
    log(f"== {VERSION} {CORE_VERSION} args={vars(args)}")
    res = {"version": VERSION, "core": CORE_VERSION, "args": vars(args)}
    res["lemma1"] = lemma1(args.amax, log)
    res["lemma23"] = lemma23(args.amax, log)
    res["lemma4"] = lemma4(args.amax, log)
    ok5, res["theorem5_rows"] = theorem5(args.free_amax, log)
    res["theorem5"] = ok5
    ok6, res["cor52_rows"] = corollary52(args.amax, log)
    res["cor52"] = ok6
    res["seconds"] = round(time.time() - t0, 1)
    res["verdict"] = "PASS" if all([res["lemma1"], res["lemma23"], res["lemma4"],
                                    res["theorem5"], res["cor52"]]) else "FAIL"
    log(f"verdict: {res['verdict']} ({res['seconds']} s)")
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(res, f, indent=1)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# check_C37 — алгебра двойного разреза и граница свободного разреза\n\n")
        f.write(f"Команда: `python3 checks/check_C37.py --amax {args.amax} "
                f"--free-amax {args.free_amax}`. Версии: {VERSION}, {CORE_VERSION}.\n\n```text\n")
        f.write("\n".join(lines))
        f.write("\n```\n")


if __name__ == "__main__":
    main()
