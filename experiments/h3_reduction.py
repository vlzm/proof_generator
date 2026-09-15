"""H3 (PLAN §4.3, item 3): is the one-element reduction compatible with the
budget n - 1?

The induction that would give C15 is D_n <= D_{n-1} + (n - 1) (then
D_4 = 6 = B_4 and B_n - B_{n-1} = n - 1 give D_n <= B_n with zero slack).
A scheme proving it must, from any state of the n-problem, reach in at most
n - 1 moves a state that is "the (n-1)-problem plus one parked element", and
then simulate a shortest word of the smaller problem move for move.

This module tests the arithmetic preconditions of any such scheme, using the
certified tables (C2v).  A reduction is `reduce(p, k)`: delete one array entry
and renormalise the values.  Two statistics per state p:

  (1) excess(p)  = d_n(p) - min over k of d_{n-1}(reduce(p,k)).
      If excess > n - 1, then that (cheapest) reduction is provably incompatible
      with a 1:1 simulation: the scheme would output a word shorter than d_n(p).
  (2) slack(p)   = max over k of d_{n-1}(reduce(p,k)) - (d_n(p) - (n - 1)).
      A necessary condition for a scheme with zero preparation moves: some
      reduction must be hard enough, slack(p) >= 0.  With m <= n - 1
      preparation moves the condition becomes
      max over (u, k), |u| <= m, of d_{n-1}(reduce(apply(p,u), k)) >= d_n(p) - m;
      it is checked exhaustively over all words u for small n.

State conventions: PROBLEM §2-§3 (array read from the head, head at indices
(0,1)).  Reduction `reduce(p, k)`: delete the entry at array index k and shift
the values above p[k] down by one; the head stays at index 0 of the new array
(for k = 0 the array starts at the old index 1, i.e. the head moves with it).

Usage: python3 experiments/h3_reduction.py [--nmax 8] [--sample-from 9] [--sample 200000]
Output: data/runs/h3_reduction/report.md, report.json.
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
sys.path.insert(0, os.path.join(ROOT, "exact"))

from bfs import factorials, rank_perm  # noqa: E402
from moves import CORE_VERSION, apply_word  # noqa: E402

VERSION = "h3_reduction-1.0"
OUT = os.path.join(ROOT, "data", "runs", "h3_reduction")
TABLES = os.path.join(ROOT, "data", "tables")


def load_table(n):
    with open(os.path.join(TABLES, f"dist_n{n}.bin"), "rb") as f:
        return f.read()


def all_words(maxlen):
    """All words over {L, R, X} of length 0..maxlen, shortest first."""
    out = [""]
    cur = [""]
    for _ in range(maxlen):
        cur = [w + c for w in cur for c in "LRX"]
        out += cur
    return out


def reduce_state(p, k):
    """Delete array entry k, renormalise values; returns a tuple of length n-1."""
    v = p[k]
    return tuple(x - 1 if x > v else x for i, x in enumerate(p) if i != k)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8)
    ap.add_argument("--sample-from", type=int, default=9)
    ap.add_argument("--sample", type=int, default=200000)
    ap.add_argument("--prep-max", type=int, default=7)
    ap.add_argument("--seed", type=int, default=20260915)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    rng = random.Random(args.seed)
    lines = []
    res = {"version": VERSION, "core": CORE_VERSION, "args": vars(args), "n": {}}

    def log(m):
        print(m)
        lines.append(m)

    log(f"== {VERSION} core={CORE_VERSION} args={vars(args)}")
    ok = True
    top = max(args.nmax, args.sample_from if args.sample else args.nmax)
    for n in range(5, top + 1):
        t0 = time.time()
        try:
            tab_n = load_table(n)
            tab_m = load_table(n - 1)
        except FileNotFoundError:
            log(f"n={n}: table missing, skipped")
            continue
        fact_n = factorials(n)
        fact_m = factorials(n - 1)
        exhaustive = n <= args.nmax
        if exhaustive:
            src = itertools.permutations(range(n))
            cover = "exhaustive"
            total = fact_n[n]
        else:
            total = args.sample
            src = (tuple(rng.sample(range(n), n)) for _ in range(total))
            cover = f"random sample of {total}"
        worst = -10 ** 9
        worst_p = None
        hist = {}
        # per-fixed-rule maxima: rule k = "always delete array index k"
        rule_max = [-10 ** 9] * n
        cnt = 0
        min_slack = 10 ** 9
        min_slack_p = None
        neg_slack = 0
        for p in src:
            d = tab_n[rank_perm(p, fact_n)]
            best = None
            hard = None
            for k in range(n):
                dm = tab_m[rank_perm(reduce_state(p, k), fact_m)]
                if best is None or dm < best:
                    best = dm
                if hard is None or dm > hard:
                    hard = dm
                if d - dm > rule_max[k]:
                    rule_max[k] = d - dm
            e = d - best
            hist[e] = hist.get(e, 0) + 1
            if e > worst:
                worst, worst_p = e, p
            sl = hard - (d - (n - 1))
            if sl < 0:
                neg_slack += 1
            if sl < min_slack:
                min_slack, min_slack_p = sl, p
            cnt += 1
        verdict = ("<= n-1" if worst <= n - 1 else
                   "> n-1: the cheapest reduction is incompatible with a 1:1 simulation")
        log(f"n={n} ({cover}, {cnt} states, tables C2v): "
            f"max [d_n(p) - min_k d_(n-1)(reduce(p,k))] = {worst} vs n-1 = {n-1} "
            f"-> {verdict}; worst p = {worst_p}; {time.time()-t0:.1f} s")
        log(f"   histogram of the excess: "
            + ", ".join(f"{e}:{hist[e]}" for e in sorted(hist)))
        log(f"   best fixed deletion rule (delete array index k): "
            + ", ".join(f"k={k}:{rule_max[k]}" for k in range(n))
            + f"  (budget {n-1})")
        log(f"   slack (zero preparation): min over states of "
            f"[max_k d_(n-1) - (d_n - (n-1))] = {min_slack} at p = {min_slack_p}; "
            f"states with negative slack: {neg_slack} of {cnt}")
        if neg_slack:
            ok = False
        res["n"][n] = {"coverage": cover, "states": cnt, "worst_excess": worst,
                       "worst_p": list(worst_p), "budget": n - 1,
                       "histogram": {str(k): v for k, v in sorted(hist.items())},
                       "fixed_rule_max": rule_max,
                       "min_slack": min_slack, "min_slack_p": list(min_slack_p),
                       "negative_slack_states": neg_slack}
        if exhaustive and n <= args.prep_max:
            # can <= n-1 preparation moves repair the negative-slack states?
            words = all_words(n - 1)
            rescued = 0
            failed = []
            for p in itertools.permutations(range(n)):
                d = tab_n[rank_perm(p, fact_n)]
                hard = max(tab_m[rank_perm(reduce_state(p, k), fact_m)]
                           for k in range(n))
                if hard - (d - (n - 1)) >= 0:
                    continue
                good = False
                for u in words:
                    q = apply_word(p, u)
                    hq = max(tab_m[rank_perm(reduce_state(q, k), fact_m)]
                             for k in range(n))
                    if hq >= d - len(u):
                        good = True
                        break
                if good:
                    rescued += 1
                else:
                    failed.append(p)
            log(f"   preparation (all words |u| <= {n-1}): rescued {rescued} of "
                f"{neg_slack} negative-slack states; unrescued: {len(failed)}"
                + (f", e.g. {failed[:3]}" if failed else ""))
            res["n"][n]["prep_rescued"] = rescued
            res["n"][n]["prep_unrescued"] = len(failed)
            res["n"][n]["prep_unrescued_examples"] = [list(x) for x in failed[:5]]
    log(f"verdict: {'no obstruction found (slack >= 0 everywhere tested)' if ok else 'OBSTRUCTION: states with negative slack exist'}")
    res["verdict"] = "no-obstruction" if ok else "obstruction"
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(res, f, indent=1)
    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# h3_reduction — индукция H3: бюджет n−1 против одноэлементной редукции\n\n")
        f.write(f"Команда: `python3 experiments/h3_reduction.py --nmax {args.nmax} "
                f"--sample-from {args.sample_from} --sample {args.sample} "
                f"--prep-max {args.prep_max} --seed {args.seed}`. "
                f"Версии: {VERSION}, {CORE_VERSION}, таблицы C2v.\n\n```text\n")
        f.write("\n".join(lines))
        f.write("\n```\n")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
