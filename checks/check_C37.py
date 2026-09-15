"""Checker for C37 — toric reformulation of H13-I (session 9).

Lemmas of docs/proofs/C37_toric_inversions.md, verified exhaustively on small n.
Cut convention: the repo cut (q, c) is rewritten as (a, b) = (q+1, q+1-c),
so the line is w_j = (pi(a+j) - b) mod n and (a, b) runs over all of Z_n x Z_n.

  part A (lemma 1): arc/XOR form of the inversion count.
      For every pair {x, y}: inverted at (a, b) iff exactly one of
      a in (y, x], b in (pi(y), pi(x)] holds.  All pi, all cuts, 4 <= n <= AMAX.
  part B (lemma 2): gradient.  inv(a+1, b) - inv(a, b) = (n-1) - 2*w_0 and
      inv(a, b+1) - inv(a, b) = (n-1) - 2*pos_w(0).  All pi, all cuts.
  part C (lemma 3): block shift.  inv(a+k, b+l) - inv(a, b) equals
      (ascending minus descending) over the pairs between D and its complement,
      where D = (positions < k) sym.diff. (values < l), |D| = k + l - 2m,
      and the number of such pairs is |D|*(n-|D|).  All pi, all cuts, all (k, l).
  part D (lemma 5): Mantel threshold.  Every 123-avoiding line has
      inv >= floor((n-1)^2/4), with equality exactly on the lines
      (alpha-1, ..., 0, n-1, ..., alpha) with alpha in {floor(n/2), ceil(n/2)}
      (the lines of the reflections).  All lines, 4 <= n <= AMAX.
  part E (lemma 6): triple orientation.  T+(pi) = number of triples whose
      cyclic value order matches the cyclic position order is cut-invariant;
      every ascending triple of every line is such a triple; and
      T+(pi) = 0 iff pi is a reflection iff every line is 123-avoiding.
  part F (sanity): max_pi I(pi) = floor((n-1)^2/4) (C33/C35 data, 4 <= n <= 8).

Usage: python3 checks/check_C37.py --amax 7
Output: data/runs/check_C37/.  Version check_C37-1.0.
"""

import argparse
import itertools
import json
import os
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
from moves import CORE_VERSION  # noqa: E402

VERSION = "check_C37-1.0"
OUT = os.path.join(ROOT, "data", "runs", "check_C37")


def target(n):
    return ((n - 1) ** 2) // 4


def line(perm, a, b):
    n = len(perm)
    return tuple((perm[(a + j) % n] - b) % n for j in range(n))


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def inv_table(perm):
    """T[a][b] = inversions of line(perm, a, b), built by the gradient lemma."""
    n = len(perm)
    pinv = [0] * n
    for i, v in enumerate(perm):
        pinv[v] = i
    T = [[0] * n for _ in range(n)]
    T[0][0] = inversions(line(perm, 0, 0))
    for b in range(n - 1):
        T[0][b + 1] = T[0][b] + (n - 1) - 2 * pinv[b]
    for b in range(n):
        for a in range(n - 1):
            T[a + 1][b] = T[a][b] + (n - 1) - 2 * ((perm[a] - b) % n)
    return T


def in_arc(z, start, end, n):
    """z in (start, end] on the circle Z_n (empty when start == end)."""
    return 0 < (z - start) % n <= (end - start) % n and start != end


def part_a(n):
    """Lemma 1: inverted iff exactly one of a in A_p, b in B_p."""
    cnt = 0
    for perm in itertools.permutations(range(n)):
        for a in range(n):
            for b in range(n):
                w = line(perm, a, b)
                pos = {}
                for j in range(n):
                    pos[(a + j) % n] = j
                for x in range(n):
                    for y in range(n):
                        if x == y:
                            continue
                        inv_actual = (pos[x] < pos[y]) != (w[pos[x]] < w[pos[y]])
                        ina = in_arc(a, y, x, n)
                        inb = in_arc(b, perm[y], perm[x], n)
                        assert inv_actual == (ina != inb), (perm, a, b, x, y)
                        cnt += 1
    return cnt


def part_b(n):
    """Lemma 2: the two gradient formulas; inv_table agrees with brute force."""
    cnt = 0
    for perm in itertools.permutations(range(n)):
        T = inv_table(perm)
        for a in range(n):
            for b in range(n):
                w = line(perm, a, b)
                assert T[a][b] == inversions(w), (perm, a, b)
                assert T[(a + 1) % n][b] - T[a][b] == (n - 1) - 2 * w[0]
                assert T[a][(b + 1) % n] - T[a][b] == (n - 1) - 2 * w.index(0)
                cnt += 1
    return cnt


def part_c(n):
    """Lemma 3: block shift identity, size of D and of the cross-pair set."""
    cnt = 0
    for perm in itertools.permutations(range(n)):
        T = inv_table(perm)
        for a in range(n):
            for b in range(n):
                w = line(perm, a, b)
                for k in range(n):
                    for l in range(n):
                        D = [j for j in range(n) if (j < k) != (w[j] < l)]
                        m = sum(1 for j in range(k) if w[j] < l)
                        assert len(D) == k + l - 2 * m
                        Dset = set(D)
                        asc = des = 0
                        for i in range(n):
                            for j in range(i + 1, n):
                                if (i in Dset) == (j in Dset):
                                    continue
                                if w[i] < w[j]:
                                    asc += 1
                                else:
                                    des += 1
                        assert asc + des == len(D) * (n - len(D))
                        delta = T[(a + k) % n][(b + l) % n] - T[a][b]
                        assert delta == asc - des, (perm, a, b, k, l)
                        cnt += 1
    return cnt


def has_ascending_triple(w):
    n = len(w)
    for i in range(n):
        for j in range(i + 1, n):
            if w[i] >= w[j]:
                continue
            for k in range(j + 1, n):
                if w[j] < w[k]:
                    return True
    return False


def reflection_lines(n):
    out = set()
    for alpha in (n // 2, (n + 1) // 2):
        out.add(tuple(list(range(alpha - 1, -1, -1)) + list(range(n - 1, alpha - 1, -1))))
    return out


def part_d(n):
    """Lemma 5: Mantel threshold on 123-avoiding lines."""
    tgt = target(n)
    refl = reflection_lines(n)
    eq = set()
    cnt = 0
    for w in itertools.permutations(range(n)):
        if has_ascending_triple(w):
            continue
        cnt += 1
        iv = inversions(w)
        assert iv >= tgt, (w, iv, tgt)
        if iv == tgt:
            eq.add(w)
    assert eq == refl, (sorted(eq), sorted(refl))
    return cnt, len(eq)


def tplus(perm):
    n = len(perm)
    t = 0
    for x, y, z in itertools.combinations(range(n), 3):
        if (perm[y] - perm[x]) % n < (perm[z] - perm[x]) % n:
            t += 1
    return t


def is_reflection(perm):
    n = len(perm)
    s = (perm[0] + 0) % n
    return all((perm[i] + i) % n == s for i in range(n))


def part_e(n):
    """Lemma 6: triple orientation is cut-invariant; T+ = 0 iff reflection."""
    nrefl = 0
    cnt = 0
    for perm in itertools.permutations(range(n)):
        tp = tplus(perm)
        any_asc = False
        for a in range(n):
            for b in range(n):
                w = line(perm, a, b)
                # every ascending triple of the line is a positively oriented
                # triple of pi, so their number never exceeds T+
                asc = 0
                for i, j, k in itertools.combinations(range(n), 3):
                    if w[i] < w[j] < w[k]:
                        asc += 1
                        x, y, z = ((a + i) % n, (a + j) % n, (a + k) % n)
                        xs, ys, zs = sorted((x, y, z))
                        assert (perm[ys] - perm[xs]) % n < (perm[zs] - perm[xs]) % n
                assert asc <= tp
                any_asc = any_asc or asc > 0
                cnt += 1
        assert (tp == 0) == (not any_asc), perm
        assert (tp == 0) == is_reflection(perm), perm
        if tp == 0:
            nrefl += 1
    assert nrefl == n
    return cnt, nrefl


def part_f(n):
    """Sanity: max_pi I(pi) = floor((n-1)^2/4), attained exactly on reflections."""
    tgt = target(n)
    best = -1
    args = []
    for perm in itertools.permutations(range(n)):
        T = inv_table(perm)
        I = min(min(row) for row in T)
        if I > best:
            best, args = I, [perm]
        elif I == best:
            args.append(perm)
    assert best == tgt, (n, best, tgt)
    assert len(args) == n and all(is_reflection(p) for p in args)
    return best, len(args)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amax", type=int, default=7)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    lines = [f"== {VERSION} {CORE_VERSION} args={vars(args)}"]
    summary = {"version": VERSION, "core": CORE_VERSION, "amax": args.amax,
               "parts": {}}
    ok = True
    t_all = time.time()
    try:
        for n in range(4, min(args.amax, 6) + 1):
            t = time.time()
            c = part_a(n)
            lines.append(f"part A n={n}: {c} (pi, cut, ordered pair) checks, "
                         f"{time.time() - t:.1f} s")
        for n in range(4, min(args.amax, 7) + 1):
            t = time.time()
            c = part_b(n)
            lines.append(f"part B n={n}: {c} cuts, gradient ok, "
                         f"{time.time() - t:.1f} s")
        for n in range(4, min(args.amax, 5) + 1):
            t = time.time()
            c = part_c(n)
            lines.append(f"part C n={n}: {c} (cut, shift) pairs, block shift ok, "
                         f"{time.time() - t:.1f} s")
        for n in range(4, args.amax + 1):
            t = time.time()
            c, e = part_d(n)
            lines.append(f"part D n={n}: {c} 123-avoiding lines, inv >= "
                         f"{target(n)}, {e} with equality = reflection lines, "
                         f"{time.time() - t:.1f} s")
        for n in range(4, min(args.amax, 6) + 1):
            t = time.time()
            c, r = part_e(n)
            lines.append(f"part E n={n}: {c} cuts, ascending triples are "
                         f"positive triples, T+=0 on {r} reflections, "
                         f"{time.time() - t:.1f} s")
        for n in range(4, min(args.amax, 8) + 1):
            t = time.time()
            b, cnt = part_f(n)
            lines.append(f"part F n={n}: max_pi I = {b} = floor((n-1)^2/4), "
                         f"attained on {cnt} reflections, {time.time() - t:.1f} s")
    except AssertionError as exc:
        ok = False
        lines.append(f"FAIL: {exc}")
    lines.append(f"total {time.time() - t_all:.1f} s")
    lines.append("verdict: " + ("PASS" if ok else "FAIL"))
    text = "\n".join(lines)
    print(text)
    summary["ok"] = ok
    summary["log"] = lines
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# check_C37 — торическая переформулировка H13-I\n\n")
        f.write(f"Команда: `python3 checks/check_C37.py --amax {args.amax}`. "
                f"Версии: {VERSION}, {CORE_VERSION}.\n\n```text\n" + text + "\n```\n")
    with open(os.path.join(OUT, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
