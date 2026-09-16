"""Algebra of the double cut (H13-I): exact identities, the free-cut bound, and
the local-maximum criterion.

Setting (docs/notes/h13_line_model.md).  For a permutation pi of Z_n a double
cut (q, c) produces the line w_j = (pi(q + 1 + j) - (q + 1 - c)) mod n.  We use
the equivalent parameters m = (q + 1) mod n (the first position of the line) and
k = (q + 1 - c) mod n (the value that becomes 0); (q, c) <-> (m, k) is a
bijection of Z_n x Z_n, so I(pi) = min_{m,k} inv(w(m, k)) is unchanged.

Objects (all for the point set {(i, pi(i))}):
  Inv(pi)  set of inverted pairs;   inv(pi) = |Inv(pi)|
  A_m      points with position < m ("position prefix")
  B_k      points with value < k    ("value prefix")
  S(m, k) = A_m xor B_k            (points of the two "anti-diagonal" quadrants)
  delta(S) set of pairs separated by S (the cut of K_n induced by S)
  e_inv(S) number of inverted pairs separated by S

Part 1 (parametrisation)  the (q, c) and (m, k) families of lines coincide.
Part 2 (identity L1)      inv(w(m, k)) = inv(pi) + s(n - s) - 2 e_inv(S),
                          s = |S|; equivalently inv(w) = |Inv(pi) xor delta(S)|.
Part 3 (pair count L2)    a pair at clockwise distances (d_p, d_v) is concordant
                          for exactly (n - d_p)(n - d_v) + d_p d_v of the n^2
                          cuts; summing gives Sigma(pi) = sum over cuts of the
                          number of concordant pairs.
Part 4 (averaging)        mean number of concordant pairs = Sigma(pi) / n^2;
                          closed forms for id and sigma_n; how many pi are
                          settled by averaging alone.
Part 5 (free cut L3)      max over ALL subsets S of [n] of the number of
                          concordant pairs is >= floor(n^2/4) (proof:
                          docs/proofs/C37_double_cut_algebra.md), and how far
                          the toric family is from that free optimum.
Part 6 (criterion H14)    is there a toric cut whose line has every element
                          inverted with at most floor((n-1)/2) others (plus, for
                          odd n, the pairwise condition on the extremal set)?
                          Such a cut would prove H13-I by L3's argument.
Part 7 (H13-I)            independent re-check of I(pi) <= floor((n-1)^2/4).

Usage: python3 experiments/toric_cut_algebra.py [--nmax 8] [--free-nmax 7]
Output: data/runs/h13i_verdict/report.json, report.md.
Version toric_cut-1.0.
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

VERSION = "toric_cut-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h13i_verdict")


# ---------------------------------------------------------------- basic tools

def pc(x):
    return bin(x).count("1")


def line(perm, m, k):
    """The line of the double cut (m, k): first position m, value k becomes 0."""
    n = len(perm)
    return [(perm[(m + j) % n] - k) % n for j in range(n)]


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def inv_masks(perm):
    """inv_masks[i] = bitmask of j != i forming an inverted pair with i."""
    n = len(perm)
    im = [0] * n
    for i in range(n):
        for j in range(n):
            if i != j and (i < j) != (perm[i] < perm[j]):
                im[i] |= 1 << j
    return im


def cut_set(perm, m, k):
    n = len(perm)
    S = 0
    for i in range(n):
        if (i < m) != (perm[i] < k):
            S |= 1 << i
    return S


def inv_from_set(im, inv0, S, n):
    """inv(pi) + s(n-s) - 2 e_inv(S)  =  |Inv(pi) xor delta(S)|."""
    full = (1 << n) - 1
    comp = full ^ S
    s = pc(S)
    e = 0
    for i in range(n):
        if S >> i & 1:
            e += pc(im[i] & comp)
    return inv0 + s * (n - s) - 2 * e


def per_element(im, S, n):
    """For the line given by S: number of elements inverted with each element."""
    full = (1 << n) - 1
    comp = full ^ S
    out = []
    for p in range(n):
        cross = (comp if (S >> p & 1) else S) & ~(1 << p)
        out.append(pc(im[p] ^ cross))
    return out


def toric_sets(perm):
    n = len(perm)
    return [(m, k, cut_set(perm, m, k)) for m in range(n) for k in range(n)]


# ------------------------------------------------------------------ part 1, 2

def part_param(nmax, log):
    ok = True
    for n in range(4, min(nmax, 6) + 1):
        cnt = 0
        for perm in itertools.permutations(range(n)):
            a = set()
            for q in range(n):
                for c in range(n):
                    sh = q + 1 - c
                    a.add(tuple((perm[(q + 1 + j) % n] - sh) % n for j in range(n)))
            b = set(tuple(line(perm, m, k)) for m in range(n) for k in range(n))
            if a != b:
                ok = False
            cnt += 1
        log(f"part 1 n={n}: {cnt} perms, (q,c) and (m,k) families of lines coincide -> {'ok' if ok else 'FAIL'}")
    return ok


def part_identity(nmax, log, samples=200, big=(8, 12, 20, 33, 50)):
    ok = True
    rows = []
    for n in range(4, nmax + 1):
        t0 = time.time()
        cnt = 0
        for perm in itertools.permutations(range(n)):
            im = inv_masks(perm)
            inv0 = sum(map(pc, im)) // 2
            for m in range(n):
                for k in range(n):
                    S = cut_set(perm, m, k)
                    if inv_from_set(im, inv0, S, n) != inversions(line(perm, m, k)):
                        ok = False
            cnt += 1
        rows.append({"n": n, "perms": cnt, "cuts": n * n, "time": round(time.time() - t0, 1)})
        log(f"part 2 n={n}: {cnt} perms x {n * n} cuts, identity L1 holds -> {'ok' if ok else 'FAIL'}"
            f" ({rows[-1]['time']} s)")
    rnd = random.Random(20260916)
    for n in big:
        for _ in range(samples):
            perm = list(range(n))
            rnd.shuffle(perm)
            im = inv_masks(perm)
            inv0 = sum(map(pc, im)) // 2
            m = rnd.randrange(n)
            k = rnd.randrange(n)
            S = cut_set(perm, m, k)
            if inv_from_set(im, inv0, S, n) != inversions(line(perm, m, k)):
                ok = False
        log(f"part 2 n={n}: {samples} random (pi, m, k), identity L1 holds -> {'ok' if ok else 'FAIL'}")
    return ok, rows


def part_pairs(nmax, log):
    """Pair concordance count over the n^2 cuts, and Sigma(pi)."""
    ok = True
    for n in range(4, min(nmax, 6) + 1):
        for perm in itertools.permutations(range(n)):
            total = 0
            for i in range(n):
                for j in range(i + 1, n):
                    d_p = (j - i) % n
                    d_v = (perm[j] - perm[i]) % n
                    pred = (n - d_p) * (n - d_v) + d_p * d_v
                    act = 0
                    for m in range(n):
                        for k in range(n):
                            pi_ = (i - m) % n
                            pj = (j - m) % n
                            vi = (perm[i] - k) % n
                            vj = (perm[j] - k) % n
                            if (pi_ < pj) == (vi < vj):
                                act += 1
                    if pred != act:
                        ok = False
                    total += pred
            sigma = sum(n * (n - 1) // 2 - inversions(line(perm, m, k))
                        for m in range(n) for k in range(n))
            if sigma != total:
                ok = False
        log(f"part 3 n={n}: pair formula (n-d_p)(n-d_v)+d_p d_v and Sigma(pi) -> {'ok' if ok else 'FAIL'}")
    return ok


def sigma_pi(perm):
    n = len(perm)
    tot = 0
    for i in range(n):
        for j in range(i + 1, n):
            d_p = (j - i) % n
            d_v = (perm[j] - perm[i]) % n
            tot += (n - d_p) * (n - d_v) + d_p * d_v
    return tot


def part_averaging(nmax, log):
    """Mean number of concordant pairs; closed forms; how many pi averaging settles."""
    rows = []
    ok = True
    for n in range(4, nmax + 1):
        tgt = (n * n) // 4
        idp = tuple(range(n))
        sg = tuple((1 - i) % n for i in range(n))
        rv = tuple(n - 1 - i for i in range(n))
        mean_id = sigma_pi(idp) / n ** 2
        mean_sg = sigma_pi(sg) / n ** 2
        mean_rv = sigma_pi(rv) / n ** 2
        if abs(mean_id - (n - 1) * (2 * n - 1) / 6) > 1e-9:
            ok = False
        if abs(mean_sg - (n * n - 1) / 6) > 1e-9 or abs(mean_rv - (n * n - 1) / 6) > 1e-9:
            ok = False
        good = 0
        total = 0
        for perm in itertools.permutations(range(n)):
            total += 1
            if sigma_pi(perm) >= tgt * n ** 2:
                good += 1
        rows.append({"n": n, "target": tgt, "mean_id": mean_id, "mean_sigma_n": mean_sg,
                     "settled_by_averaging": good, "perms": total})
        log(f"part 4 n={n}: target floor(n^2/4)={tgt}; mean(id)={mean_id:.3f} "
            f"(= (n-1)(2n-1)/6), mean(sigma_n)={mean_sg:.3f} (= (n^2-1)/6); "
            f"averaging settles {good}/{total} perms")
    return ok, rows


def best_free(im, inv0, n):
    """min over ALL subsets S of [n] of |Inv(pi) xor delta(S)| (S and its complement agree)."""
    best = None
    arg = 0
    for S in range(1 << (n - 1)):
        v = inv_from_set(im, inv0, S, n)
        if best is None or v < best:
            best = v
            arg = S
    return best, arg


def local_search_free(im, inv0, n, S=0):
    """Local search over all subsets (1- and 2-flips) -- the algorithm behind L3.

    Flipping p changes the inversion count by r_p = n - 1 - 2 inv_p; flipping
    {p, q} changes it by r_p + r_q - 2 N_pq, N_pq = +1 iff the pair is concordant
    in the current line.  Each accepted move strictly decreases the count, so the
    search stops after at most B_n steps.
    """
    cur = inv_from_set(im, inv0, S, n)
    while True:
        pe = per_element(im, S, n)
        r = [n - 1 - 2 * v for v in pe]
        p = min(range(n), key=lambda i: r[i])
        if r[p] < 0:
            S ^= 1 << p
            cur += r[p]
            continue
        full = (1 << n) - 1
        comp = full ^ S
        moved = False
        for p in range(n):
            linv_p = im[p] ^ ((comp if (S >> p & 1) else S) & ~(1 << p))
            for q in range(p + 1, n):
                npq = -1 if (linv_p >> q & 1) else 1
                d = r[p] + r[q] - 2 * npq
                if d < 0:
                    S ^= (1 << p) | (1 << q)
                    cur += d
                    moved = True
                    break
            if moved:
                break
        if not moved:
            return cur, S


def part_free(nmax, free_nmax, log):
    """L3 (free cut) exhaustively, and the gap between the free and the toric optimum."""
    ok = True
    rows = []
    for n in range(4, free_nmax + 1):
        tgt = ((n - 1) ** 2) // 4
        worst_free = -1
        worst_toric = -1
        gap_hist = {}
        gap_examples = []
        affine_gap = 0
        t0 = time.time()
        aff = set(tuple((a * i + b) % n for i in range(n))
                  for a in range(1, n) if gcd(a, n) == 1 for b in range(n))
        for perm in itertools.permutations(range(n)):
            im = inv_masks(perm)
            inv0 = sum(map(pc, im)) // 2
            bf, _ = best_free(im, inv0, n)
            bt = min(inv_from_set(im, inv0, S, n) for (_, _, S) in toric_sets(perm))
            if bf > tgt:
                ok = False
            worst_free = max(worst_free, bf)
            worst_toric = max(worst_toric, bt)
            gap_hist[bt - bf] = gap_hist.get(bt - bf, 0) + 1
            if bt > bf:
                if perm in aff:
                    affine_gap += 1
                if len(gap_examples) < 12:
                    gap_examples.append({"pi": list(perm), "I": bt, "free": bf,
                                         "affine": perm in aff})
        rows.append({"n": n, "target": tgt, "max_min_free": worst_free,
                     "max_min_toric": worst_toric,
                     "gap_histogram": {str(k): v for k, v in sorted(gap_hist.items())},
                     "positive_gap_affine": affine_gap,
                     "positive_gap_examples": gap_examples,
                     "time": round(time.time() - t0, 1)})
        npos = sum(v for k, v in gap_hist.items() if k > 0)
        log(f"part 5 n={n}: max_pi min_free = {worst_free} (target floor((n-1)^2/4) = {tgt}), "
            f"max_pi I(pi) = {worst_toric}; I - min_free histogram "
            f"{dict(sorted(gap_hist.items()))}; of the {npos} perms with a gap "
            f"{affine_gap} are affine -> {'ok' if ok else 'FAIL'}")
    # large n: the local search of the proof
    rnd = random.Random(20260916)
    for n in (12, 20, 33, 50, 100, 101):
        tgt = ((n - 1) ** 2) // 4
        worst = -1
        fams = families(n)
        for perm in fams + [random_perm(n, rnd) for _ in range(20)]:
            im = inv_masks(perm)
            inv0 = sum(map(pc, im)) // 2
            v, _ = local_search_free(im, inv0, n)
            worst = max(worst, v)
            if v > tgt:
                ok = False
        log(f"part 5 n={n}: local search over free subsets on {len(fams) + 20} inputs, "
            f"max = {worst} <= {tgt} -> {'ok' if ok else 'FAIL'}")
    return ok, rows


def part_h14(nmax, log):
    """Criterion H14: a toric cut whose line is a 1- and 2-flip local optimum."""
    rows = []
    for n in range(4, nmax + 1):
        h = (n - 1) // 2
        fails = []
        t0 = time.time()
        for perm in itertools.permutations(range(n)):
            im = inv_masks(perm)
            ok_i = False
            ok_both = False
            for (_, _, S) in toric_sets(perm):
                pe = per_element(im, S, n)
                if max(pe) <= h:
                    ok_i = True
                    if n % 2 == 0:
                        ok_both = True
                        break
                    full = (1 << n) - 1
                    comp = full ^ S
                    Z = [p for p in range(n) if pe[p] == h]
                    good = all((im[p] ^ (comp if (S >> p & 1) else S)) >> q & 1
                               for p in Z for q in Z if p != q)
                    if good:
                        ok_both = True
                        break
            if not ok_both:
                fails.append((list(perm), ok_i))
        aff = set(tuple((a * i + b) % n for i in range(n))
                  for a in range(1, n) if gcd(a, n) == 1 for b in range(n))
        naff = sum(1 for f in fails if tuple(f[0]) in aff)
        rows.append({"n": n, "failures": len(fails), "affine_failures": naff,
                     "examples": [f[0] for f in fails[:16]],
                     "time": round(time.time() - t0, 1)})
        log(f"part 6 n={n}: cuts failing criterion H14 for {len(fails)} perms"
            + (f" ({naff} affine); e.g. {fails[0][0]}" if fails else ""))
    return rows


def random_perm(n, rnd):
    p = list(range(n))
    rnd.shuffle(p)
    return tuple(p)


def families(n):
    """Structured inputs: reflections, affine, rotations, cycles, transposition products."""
    out = []
    for h in range(n):
        out.append(tuple((h - i) % n for i in range(n)))
    for a in range(1, n):
        if gcd(a, n) == 1:
            for b in (0, 1, n // 2):
                out.append(tuple((a * i + b) % n for i in range(n)))
    out.append(tuple(range(n)))
    out.append(tuple((i + n // 2) % n for i in range(n)))
    p = list(range(n))
    for i in range(0, n - 1, 2):
        p[i], p[i + 1] = p[i + 1], p[i]
    out.append(tuple(p))
    p = list(range(n))
    for i in range(0, n - 3, 4):
        p[i], p[i + 2] = p[i + 2], p[i]
        p[i + 1], p[i + 3] = p[i + 3], p[i + 1]
    out.append(tuple(p))
    return sorted(set(out))


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def part_h13i(nmax, log):
    """Independent re-check of H13-I on all perms and on structured families."""
    ok = True
    rows = []
    for n in range(4, nmax + 1):
        tgt = ((n - 1) ** 2) // 4
        worst = -1
        tight = 0
        t0 = time.time()
        for perm in itertools.permutations(range(n)):
            im = inv_masks(perm)
            inv0 = sum(map(pc, im)) // 2
            bt = min(inv_from_set(im, inv0, S, n) for (_, _, S) in toric_sets(perm))
            if bt > tgt:
                ok = False
            if bt == tgt:
                tight += 1
            worst = max(worst, bt)
        rows.append({"n": n, "target": tgt, "max_I": worst, "tight_count": tight,
                     "time": round(time.time() - t0, 1)})
        log(f"part 7 n={n}: max_pi I(pi) = {worst} (target {tgt}), attained on {tight} perms"
            f" -> {'ok' if ok else 'FAIL'}")
    rnd = random.Random(11)
    for n in (20, 50, 100, 101):
        tgt = ((n - 1) ** 2) // 4
        fams = families(n)
        if len(fams) > 40:
            fams = fams[::max(1, len(fams) // 40)]
        worst = -1
        for perm in fams + [random_perm(n, rnd) for _ in range(20)]:
            im = inv_masks(perm)
            inv0 = sum(map(pc, im)) // 2
            bt = min(inv_from_set(im, inv0, cut_set(perm, m, k), n)
                     for m in range(n) for k in range(n))
            worst = max(worst, bt)
            if bt > tgt:
                ok = False
        log(f"part 7 n={n}: {len(fams) + 20} structured/random inputs, max I = {worst} <= {tgt}"
            f" -> {'ok' if ok else 'FAIL'}")
    return ok, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--free-nmax", type=int, default=7)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    t0 = time.time()
    log(f"== {VERSION} {CORE_VERSION} args={vars(args)}")
    res = {"version": VERSION, "core": CORE_VERSION, "args": vars(args)}
    res["param_ok"] = part_param(args.nmax, log)
    ok, res["identity"] = part_identity(min(args.nmax, 7), log)
    res["identity_ok"] = ok
    res["pairs_ok"] = part_pairs(args.nmax, log)
    ok, res["averaging"] = part_averaging(min(args.nmax, 8), log)
    res["averaging_ok"] = ok
    ok, res["free"] = part_free(args.nmax, args.free_nmax, log)
    res["free_ok"] = ok
    res["h14"] = part_h14(args.nmax, log)
    ok, res["h13i"] = part_h13i(args.nmax, log)
    res["h13i_ok"] = ok
    res["seconds"] = round(time.time() - t0, 1)
    verdict = all([res["param_ok"], res["identity_ok"], res["pairs_ok"],
                   res["averaging_ok"], res["free_ok"], res["h13i_ok"]])
    res["verdict"] = "PASS" if verdict else "FAIL"
    log(f"verdict: {res['verdict']} ({res['seconds']} s)")
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(res, f, indent=1)
    with open(os.path.join(OUT, "report_log.md"), "w") as f:
        f.write("# toric_cut_algebra — лог прогона\n\n")
        f.write(f"Команда: `python3 experiments/toric_cut_algebra.py --nmax {args.nmax} "
                f"--free-nmax {args.free_nmax}`. Версии: {VERSION}, {CORE_VERSION}.\n\n```text\n")
        f.write("\n".join(lines))
        f.write("\n```\n")


if __name__ == "__main__":
    main()
