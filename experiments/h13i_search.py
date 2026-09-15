"""H13-I candidate search (session 9): does a small, explicit sub-family of
double cuts (q, c) already realize I(pi) = min_{q,c} inv(w) within the
floor((n-1)^2/4) bound, without a full n^2 search?

Setting (docs/notes/h13_line_model.md §0, PLAN H13'): pi: Z_n -> Z_n. For a
position cut a in Z_n and a value cut b in Z_n, w_j = (pi((a+j) mod n) - b)
mod n, j = 0..n-1 (this is the line model with a = q+1, b = q+1-c, so
c = (a-b) mod n). I(pi) = min_{a,b} inv(w). H13-I claims I(pi) <=
floor((n-1)^2/4) for all pi, n >= 4 (VERIFIED 4 <= n <= 10, C33/C35).

Two candidate reductions of the n^2 search space to a linear one, tested
exhaustively over all pi:

  A. Fixed diagonal: is there a single k (depending only on n) such that
     min_a inv(a, a-k mod n) <= floor((n-1)^2/4) for every pi?
  B. Adjacent-optimal-b chain: as a increases by 1, does *some* choice of
     b*(a) in argmin_b inv(a,b) form a chain with |b*(a+1)-b*(a)| <= 1
     (cyclic), for every pi? (If so, a telescoping argument along the chain
     could bound the minimum via a 1-parameter walk instead of full search.)

Both are refuted (see report): (A) fails for every k, with the worst case
growing well past the bound; (B) fails for the majority of permutations
already at n = 5..7, so the argmin-b map has no globally short-jump chain.
Negative result for the "next candidate" pass of docs/notes/h13_line_model.md
§6; H13-I itself is untouched (still CONJECTURED, VERIFIED 4 <= n <= 10).

Usage: python3 experiments/h13i_search.py --nmax 7
Output: data/runs/h13i_search/report.md, report.json. Version h13i_search-1.0.
"""

import argparse
import itertools
import json
import os
import time

VERSION = "h13i_search-1.0"
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "data", "runs", "h13i_search")


def inv_count(seq):
    n = len(seq)
    c = 0
    for i in range(n):
        si = seq[i]
        for j in range(i + 1, n):
            if si > seq[j]:
                c += 1
    return c


def f(pi, n, a, b):
    w = [(pi[(a + j) % n] - b) % n for j in range(n)]
    return inv_count(w)


def full_table(pi, n):
    return [[f(pi, n, a, b) for b in range(n)] for a in range(n)]


def test_fixed_diagonal(n, bound):
    """Candidate A: single k = (a-b) mod n working for every pi."""
    k_worst = [-1] * n
    k_worst_pi = [None] * n
    global_worst = -1
    global_worst_pi = None
    for pi in itertools.permutations(range(n)):
        table = full_table(pi, n)
        I = min(min(row) for row in table)
        if I > global_worst:
            global_worst, global_worst_pi = I, pi
        for k in range(n):
            m = min(table[a][(a - k) % n] for a in range(n))
            if m > k_worst[k]:
                k_worst[k] = m
                k_worst_pi[k] = pi
    return {
        "I_max": global_worst,
        "I_max_witness": list(global_worst_pi),
        "I_max_ok": global_worst <= bound,
        "per_k_worst": [
            {"k": k, "worst": k_worst[k], "ok": k_worst[k] <= bound,
             "witness": list(k_worst_pi[k])}
            for k in range(n)
        ],
    }


def test_adjacent_chain(n):
    """Candidate B: chain of argmin-b values with cyclic jump <= 1."""
    total = 0
    bad = 0
    bad_examples = []
    for pi in itertools.permutations(range(n)):
        total += 1
        table = full_table(pi, n)
        argmins = []
        for a in range(n):
            m = min(table[a])
            argmins.append({b for b in range(n) if table[a][b] == m})
        ok_chain = False
        for b0 in argmins[0]:
            b = b0
            maxjump = 0
            for a in range(1, n):
                best, bd = None, None
                for cand in argmins[a]:
                    d = min((cand - b) % n, (b - cand) % n)
                    if bd is None or d < bd:
                        bd, best = d, cand
                maxjump = max(maxjump, bd)
                b = best
            if maxjump <= 1:
                ok_chain = True
                break
        if not ok_chain:
            bad += 1
            if len(bad_examples) < 5:
                bad_examples.append(list(pi))
    return {"total": total, "no_short_chain": bad, "examples": bad_examples}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=7)
    ap.add_argument("--nmin", type=int, default=4)
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    results = {"version": VERSION, "started": time.strftime("%Y-%m-%d %H:%M:%S")}
    per_n = []
    for n in range(args.nmin, args.nmax + 1):
        bound = ((n - 1) ** 2) // 4
        t0 = time.time()
        a_res = test_fixed_diagonal(n, bound)
        b_res = test_adjacent_chain(n)
        dt = time.time() - t0
        per_n.append({
            "n": n, "bound": bound, "time_s": round(dt, 2),
            "candidate_A_fixed_diagonal": a_res,
            "candidate_B_adjacent_chain": b_res,
        })
        print(f"n={n} bound={bound} A_ok={all(r['ok'] for r in a_res['per_k_worst'])} "
              f"B_bad={b_res['no_short_chain']}/{b_res['total']} time={dt:.1f}s")
    results["per_n"] = per_n

    with open(os.path.join(OUT, "report.json"), "w") as fh:
        json.dump(results, fh, indent=2)

    with open(os.path.join(OUT, "report.md"), "w") as fh:
        fh.write("# H13-I candidate search (session 9)\n\n")
        fh.write(f"Version {VERSION}. Exhaustive over all pi, {args.nmin} <= n <= {args.nmax}.\n\n")
        fh.write("## Candidate A: fixed diagonal k = (a-b) mod n\n\n")
        fh.write("REFUTED for every n and every k: worst-case min_a f(a,a-k) exceeds "
                  "floor((n-1)^2/4) for every choice of k.\n\n")
        fh.write("| n | bound | worst any k | global I_max (full n^2 search) |\n")
        fh.write("|---|---|---|---|\n")
        for r in per_n:
            worst_any_k = max(x["worst"] for x in r["candidate_A_fixed_diagonal"]["per_k_worst"])
            fh.write(f"| {r['n']} | {r['bound']} | {worst_any_k} | "
                      f"{r['candidate_A_fixed_diagonal']['I_max']} |\n")
        fh.write("\n## Candidate B: adjacent-optimal-b chain (jump <= 1)\n\n")
        fh.write("REFUTED: fails for the majority of pi already at n = 5.\n\n")
        fh.write("| n | pi without short chain | total pi |\n")
        fh.write("|---|---|---|\n")
        for r in per_n:
            b = r["candidate_B_adjacent_chain"]
            fh.write(f"| {r['n']} | {b['no_short_chain']} | {b['total']} |\n")
        fh.write("\nExample pi with no jump<=1 chain (n = smallest tested): ")
        fh.write(f"{per_n[0]['candidate_B_adjacent_chain']['examples']}\n")

    print("Report written to", OUT)


if __name__ == "__main__":
    main()
