"""H13-I: which restricted cut families and which inductions could work (session 9).

H13-I: for every pi there is a double cut with inv(w) <= floor((n-1)^2/4).
Cut convention: (a, b) = (q+1, q+1-c), line w_j = (pi(a+j) - b) mod n.

Part 1 — anchored families.  The only cut families of size n that are invariant
under the torus action are the *anchored* ones: for fixed (p0, v0) take the n
cuts that put some element x at line position p0 with line value v0, i.e.
(a, b) = (x - p0, pi(x) - v0).  For (p0, v0) = (0, t) this is the family "the
first element of the line has value t" already refuted in
docs/notes/h13_line_model.md §1.7.  We scan all n^2 anchored families.

Part 2 — deletion induction.  Deleting k elements (contracting both the position
circle and the value circle) gives pi' on n-k points; the induction step needs
I(pi) <= I(pi') + (floor((n-1)^2/4) - floor((n-k-1)^2/4)) for some deletion.
We report max_pi min_deletion (I(pi) - I(pi')) against the allowed gap.

Part 3 — how rare good cuts are: min_pi #{cuts with inv <= floor((n-1)^2/4)}.

Usage: python3 experiments/h13i_cut_families.py --nmax 7 [--nmax-anchor 8]
Output: data/runs/h13i_cut_families/.  Version h13i_cut_families-1.0.
"""

import argparse
import itertools
import json
import os
import resource
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
from moves import CORE_VERSION  # noqa: E402

VERSION = "h13i_cut_families-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h13i_cut_families")

AFFINE_WITNESSES = [(0, 5, 2, 7, 4, 1, 6, 3), (0, 3, 6, 1, 4, 7, 2, 5),
                    (0, 2, 4, 6, 8, 1, 3, 5, 7), (0, 5, 1, 6, 2, 7, 3, 8, 4)]


def target(n):
    return ((n - 1) ** 2) // 4


def inversions(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def inv_table(perm):
    """T[a][b] = inv of the line at cut (a, b) (gradient recurrences, C37 §2)."""
    n = len(perm)
    pinv = [0] * n
    for i, v in enumerate(perm):
        pinv[v] = i
    T = [[0] * n for _ in range(n)]
    T[0][0] = inversions(tuple(perm))
    for b in range(n - 1):
        T[0][b + 1] = T[0][b] + (n - 1) - 2 * pinv[b]
    for b in range(n):
        for a in range(n - 1):
            T[a + 1][b] = T[a][b] + (n - 1) - 2 * ((perm[a] - b) % n)
    return T


def I_of(perm):
    return min(min(row) for row in inv_table(perm))


def delete(perm, xs):
    """Delete the elements xs: contract the position circle and the value circle."""
    n = len(perm)
    xs = set(xs)
    rem = [x for x in range(n) if x not in xs]
    vals = sorted(perm[x] for x in rem)
    vrank = {v: i for i, v in enumerate(vals)}
    return tuple(vrank[perm[x]] for x in rem)


def scan_anchor(n):
    """For every (p0, v0): max over pi of min over the anchored family of inv."""
    worst = [[(-1, None) for _ in range(n)] for _ in range(n)]
    for perm in itertools.permutations(range(n)):
        T = inv_table(perm)
        for p0 in range(n):
            for v0 in range(n):
                mm = min(T[(x - p0) % n][(perm[x] - v0) % n] for x in range(n))
                if mm > worst[p0][v0][0]:
                    worst[p0][v0] = (mm, perm)
    good = [(p0, v0) for p0 in range(n) for v0 in range(n)
            if worst[p0][v0][0] <= target(n)]
    best = min(((worst[p0][v0][0], (p0, v0), worst[p0][v0][1])
                for p0 in range(n) for v0 in range(n)), key=lambda t: t[0])
    return good, best


def scan_induction(n, k):
    """max over pi of min over k-deletions of (I(pi) - I(pi'))."""
    allowed = target(n) - target(n - k)
    worst = (-10 ** 9, None)
    for perm in itertools.permutations(range(n)):
        In = I_of(perm)
        need = min(In - I_of(delete(perm, xs))
                   for xs in itertools.combinations(range(n), k))
        if need > worst[0]:
            worst = (need, perm)
    return allowed, worst


def scan_rarity(n):
    """min over pi of the number of cuts with inv <= target(n)."""
    tgt = target(n)
    worst = (10 ** 9, None)
    for perm in itertools.permutations(range(n)):
        T = inv_table(perm)
        cnt = sum(1 for a in range(n) for b in range(n) if T[a][b] <= tgt)
        if cnt < worst[0]:
            worst = (cnt, perm)
    return worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=7)
    ap.add_argument("--nmax-anchor", type=int, default=7)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    log = [f"== {VERSION} {CORE_VERSION} args={vars(args)}"]
    data = {"version": VERSION, "core": CORE_VERSION, "args": vars(args)}
    t_all = time.time()

    anchor = {}
    for n in range(4, args.nmax_anchor + 1):
        t = time.time()
        good, best = scan_anchor(n)
        anchor[n] = {"good": good, "best": [best[0], list(best[1]), list(best[2])]}
        log.append(f"part 1 n={n} target={target(n)}: anchored families reaching "
                   f"the target: {len(good)} of {n * n}; best family {best[1]} "
                   f"gives {best[0]} on pi={best[2]}, {time.time() - t:.1f} s")
    data["anchor"] = anchor

    ind = {}
    for k in (1, 2, 3):
        for n in range(max(5, k + 4), args.nmax + 1):
            t = time.time()
            allowed, worst = scan_induction(n, k)
            ind[f"{n},{k}"] = {"allowed": allowed, "needed": worst[0],
                               "pi": list(worst[1])}
            log.append(f"part 2 n={n} delete {k}: allowed gap {allowed}, worst "
                       f"needed {worst[0]} on pi={worst[1]} -> "
                       f"{'ok' if worst[0] <= allowed else 'FAILS'}, "
                       f"{time.time() - t:.1f} s")
    for perm in AFFINE_WITNESSES:
        n = len(perm)
        I = I_of(perm)
        row = []
        for k in (1, 2, 3):
            best = max((I_of(delete(perm, xs)), xs)
                       for xs in itertools.combinations(range(n), k))
            row.append((k, best[0], list(best[1]), target(n) - target(n - k),
                        I - best[0]))
            log.append(f"part 2 witness {perm}: I={I}, delete {k}: best "
                       f"I(pi')={best[0]} at {best[1]}, allowed gap "
                       f"{target(n) - target(n - k)}, needed {I - best[0]} -> "
                       f"{'ok' if I - best[0] <= target(n) - target(n - k) else 'FAILS'}")
        ind[str(perm)] = {"I": I, "rows": row}
    data["induction"] = ind

    rare = {}
    for n in range(4, args.nmax_anchor + 1):
        t = time.time()
        cnt, perm = scan_rarity(n)
        rare[n] = {"min_good_cuts": cnt, "pi": list(perm)}
        log.append(f"part 3 n={n}: min over pi of #cuts with inv <= {target(n)} "
                   f"is {cnt} (of {n * n}) on pi={perm}, {time.time() - t:.1f} s")
    data["rarity"] = rare

    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    log.append(f"total {time.time() - t_all:.1f} s, peak {mem:.0f} MB")
    text = "\n".join(log)
    print(text)
    data["log"] = log
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# h13i_cut_families — ограниченные семейства разрезов и индукция\n\n")
        f.write(f"Цель: проверить, достаточно ли для H13-I какого-нибудь якорного "
                f"семейства из n разрезов, и проходит ли индукция удалением "
                f"k элементов. Покрытие: исчерпывающе по всем pi при "
                f"4 <= n <= {args.nmax_anchor} (часть 1, 3) и "
                f"5 <= n <= {args.nmax} (часть 2), плюс именованные аффинные "
                f"свидетели при n = 8, 9. Seed не используется (перебор).\n\n")
        f.write(f"Команда: `python3 experiments/h13i_cut_families.py --nmax "
                f"{args.nmax} --nmax-anchor {args.nmax_anchor}`. "
                f"Версии: {VERSION}, {CORE_VERSION}.\n\n```text\n" + text + "\n```\n")
    with open(os.path.join(OUT, "summary.json"), "w") as f:
        json.dump(data, f, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
