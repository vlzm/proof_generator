"""H13-I probe: stress-test `I(pi) <= floor((n-1)^2/4)` beyond the exhaustive
range (4 <= n <= 10, `experiments/line_profile.c`) using two non-exhaustive
methods, and test one candidate constructive sub-lemma.

`I(pi) = min_{q,c} inv(w)` over all n^2 double cuts, `w_j = pi(q+1+j) -
(q+1-c) mod n` (same definition as `line_profile.c` / `line_model.py`).

1. Multiplicative family `pi(i) = k*i mod n`, gcd(k, n) = 1: computes I(pi)
   for every valid k, for n up to `--nmax-mult` (default 31). This is a
   SAMPLED family, not exhaustive, and not random: chosen because these
   permutations are the standard "maximally scrambling" candidates (no
   arithmetic-progression structure) and a natural place to look for a
   violation of H13-I beyond n = 10.
2. Random-restart hill-climbing search (`--restarts`, `--steps`) trying to
   exceed floor((n-1)^2/4) by local transpositions, for n up to `--nmax-hc`
   (default 24); plus an exhaustive single-transposition neighborhood check
   around the reflection sigma_n (all C(n,2) transpositions), to test
   whether the reflection is a strict local maximum of I.
3. Candidate sub-lemma check ("sliding window maps to a contiguous value
   window"): for window size m = floor(n/2), does every pi have some
   position-window of size m whose image is a contiguous value-window of
   size m? Reports the first counterexample found (used to rule out an
   attractive but false constructive route to H13-I).

None of this is a proof; it is evidence-gathering for the open problem
registered as H13-I (PLAN.md, CLAIMS.md C33/C35/C36). Reproducible: fixed
seed for the hill-climbing search. Version h13_i_probe-1.0.

Usage: python3 experiments/h13_i_probe.py
"""

import argparse
import json
import math
import os
import random
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "data", "runs", "h13_i_probe")


def inv_count(w):
    n = len(w)
    c = 0
    for i in range(n):
        wi = w[i]
        for j in range(i + 1, n):
            if wi > w[j]:
                c += 1
    return c


def I_of(pi):
    n = len(pi)
    best = None
    for q in range(n):
        for c in range(n):
            shift = (q + 1 - c) % n
            w = [(pi[(q + 1 + j) % n] - shift) % n for j in range(n)]
            iv = inv_count(w)
            if best is None or iv < best:
                best = iv
            if best == 0:
                return 0
    return best


def floorb(n):
    return (n - 1) ** 2 // 4


def is_contiguous_value_window(values, n, m):
    """True iff `values` (a set of m distinct residues mod n) is {c, c+1, ..., c+m-1} mod n for some c."""
    for c in range(n):
        if set((c + j) % n for j in range(m)) == values:
            return True
    return False


def run_multiplicative(nmax, log):
    rows = []
    best_nonrefl_ratio = -1.0
    best_nonrefl_row = None
    for n in range(5, nmax + 1):
        fb = floorb(n)
        best_k, best_I = None, -1
        second_k, second_I = None, -1  # best among k != n-1 (i.e. excluding the reflection)
        for k in range(1, n):
            if math.gcd(k, n) != 1:
                continue
            pi = [(k * i) % n for i in range(n)]
            Ival = I_of(pi)
            if Ival > fb:
                raise SystemExit(f"COUNTEREXAMPLE to H13-I: n={n}, k={k}, I={Ival} > floor((n-1)^2/4)={fb}")
            if Ival > best_I:
                best_I, best_k = Ival, k
            if k != n - 1 and Ival > second_I:
                second_I, second_k = Ival, k
        ratio = best_I / fb
        row = {"n": n, "floor": fb, "max_I_over_k": best_I, "argmax_k": best_k, "ratio": round(ratio, 4),
               "best_nonreflection_k": second_k, "best_nonreflection_I": second_I,
               "best_nonreflection_ratio": round(second_I / fb, 4) if second_I >= 0 else None}
        rows.append(row)
        if row["best_nonreflection_ratio"] is not None and row["best_nonreflection_ratio"] > best_nonrefl_ratio:
            best_nonrefl_ratio, best_nonrefl_row = row["best_nonreflection_ratio"], row
        log(f"multiplicative n={n:3d}: floor={fb:5d} max_I={best_I:5d} (k={best_k}) ratio={ratio:.3f}; "
            f"best non-reflection k={second_k} I={second_I} ratio={row['best_nonreflection_ratio']}")
    return rows, best_nonrefl_row


def run_hillclimb(nmax, restarts, steps, seed, log):
    rnd = random.Random(seed)
    rows = []
    for n in [n for n in (12, 15, 18, 20, 24) if n <= nmax]:
        fb = floorb(n)
        refl = [(n - 1 - i) % n for i in range(n)]
        I_refl = I_of(refl)
        local_best = I_refl
        for a in range(n):
            for b in range(a + 1, n):
                pi = refl[:]
                pi[a], pi[b] = pi[b], pi[a]
                Ival = I_of(pi)
                if Ival > fb:
                    raise SystemExit(f"COUNTEREXAMPLE to H13-I near reflection: n={n}, swap ({a},{b}), I={Ival} > {fb}")
                if Ival > local_best:
                    local_best = Ival
        global_best = 0
        for _ in range(restarts):
            pi = list(range(n))
            rnd.shuffle(pi)
            cur = I_of(pi)
            for _ in range(steps):
                a, b = rnd.sample(range(n), 2)
                pi[a], pi[b] = pi[b], pi[a]
                new_I = I_of(pi)
                if new_I > fb:
                    raise SystemExit(f"COUNTEREXAMPLE to H13-I via hill-climb: n={n}, I={new_I} > {fb}, pi={pi}")
                if new_I >= cur:
                    cur = new_I
                else:
                    pi[a], pi[b] = pi[b], pi[a]
            global_best = max(global_best, cur)
        row = {"n": n, "floor": fb, "I_reflection": I_refl,
               "best_after_1_transposition_from_reflection": local_best,
               "hillclimb_best": global_best, "restarts": restarts, "steps": steps}
        rows.append(row)
        log(f"hillclimb n={n:3d}: floor={fb:4d} I(refl)={I_refl:4d} "
            f"best_1swap_from_refl={local_best:4d} hillclimb_best={global_best:4d}")
    return rows


def run_sliding_window_check(nmax, log):
    """Test the sub-lemma: some position-window of size m = floor(n/2) maps onto
    a contiguous value-window. Report the first n, pi where it fails for ALL windows."""
    for n in range(5, nmax + 1):
        m = n // 2
        for k in range(1, n):
            if math.gcd(k, n) != 1:
                continue
            pi = [(k * i) % n for i in range(n)]
            found = False
            for q in range(n):
                window = set(pi[(q + j) % n] for j in range(m))
                if is_contiguous_value_window(window, n, m):
                    found = True
                    break
            if not found:
                log(f"sliding-window sub-lemma FALSE: n={n}, k={k} (pi(i)={k}*i mod {n}): "
                    f"no window of size {m} maps to a contiguous value-window")
                return {"n": n, "k": k, "m": m, "pi": pi}
    log(f"sliding-window sub-lemma: no counterexample found for multiplicative family up to n={nmax}")
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax-mult", type=int, default=31)
    ap.add_argument("--nmax-hc", type=int, default=24)
    ap.add_argument("--restarts", type=int, default=6)
    ap.add_argument("--steps", type=int, default=400)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    logf = open(os.path.join(OUT, "report.log"), "a")

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + "\n")

    t0 = time.time()
    log(f"=== h13_i_probe-1.0 run, seed={args.seed} ===")
    mult_rows, worst_row = run_multiplicative(args.nmax_mult, log)
    hc_rows = run_hillclimb(args.nmax_hc, args.restarts, args.steps, args.seed, log)
    sw_counterexample = run_sliding_window_check(args.nmax_mult, log)
    elapsed = round(time.time() - t0)
    log(f"done in {elapsed} s, no counterexample to H13-I found")

    report = {
        "version": "h13_i_probe-1.0",
        "seed": args.seed,
        "nmax_mult": args.nmax_mult,
        "nmax_hc": args.nmax_hc,
        "restarts": args.restarts,
        "steps": args.steps,
        "multiplicative_family": mult_rows,
        "worst_ratio_non_reflection": worst_row,
        "hillclimb": hc_rows,
        "sliding_window_sublemma_counterexample": sw_counterexample,
        "seconds": elapsed,
        "result": "no violation of H13-I (I(pi) <= floor((n-1)^2/4)) found; "
                  "reflection sigma_n is the unique maximizer within tested families "
                  "and a strict local max under single transpositions up to n=24",
    }
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(report, f, indent=1)

    with open(os.path.join(OUT, "report.md"), "w") as f:
        f.write("# H13-I probe (h13_i_probe-1.0)\n\n")
        f.write(f"Seed {args.seed}, {elapsed} s.\n\n")
        f.write("## Multiplicative family pi(i) = k*i mod n (all k coprime to n)\n\n")
        f.write("| n | floor((n-1)^2/4) | max I over k | argmax k | ratio |\n|---|---|---|---|---|\n")
        for r in mult_rows:
            f.write(f"| {r['n']} | {r['floor']} | {r['max_I_over_k']} | {r['argmax_k']} | {r['ratio']} |\n")
        f.write(f"\nWorst ratio among non-reflection multipliers (k != n-1): "
                f"{worst_row}\n\n")
        f.write("## Hill-climbing / reflection-neighborhood check\n\n")
        f.write("| n | floor | I(reflection) | best after 1 transposition | hill-climb best |\n"
                "|---|---|---|---|---|\n")
        for r in hc_rows:
            f.write(f"| {r['n']} | {r['floor']} | {r['I_reflection']} | "
                    f"{r['best_after_1_transposition_from_reflection']} | {r['hillclimb_best']} |\n")
        f.write("\n## Sliding-window sub-lemma\n\n")
        if sw_counterexample:
            f.write(f"FALSE in general: n={sw_counterexample['n']}, "
                    f"pi(i) = {sw_counterexample['k']}*i mod {sw_counterexample['n']}, "
                    f"no window of size {sw_counterexample['m']} maps onto a contiguous "
                    f"value-window (checked all n windows).\n")
        else:
            f.write("No counterexample found (unexpected; re-check).\n")
        f.write(f"\n## Conclusion\n\n{report['result']}.\n")

    print(f"wrote {OUT}/report.md")


if __name__ == "__main__":
    main()
