"""Independent check of C28 (docs/proofs/C28_reflections.md): the closed-form
proof that on reflections pi_h(i) = (h-i) mod n, for all n >= 4 and all
h in Z_n,

    min_c [2F_c - S_c + H(c)] <= B_n + 1     (with equality iff n = 2 mod 4, h even)
    min_c [2F_c - S_c + R_c] <= B_n + 1      (a fortiori, since R_c <= H(c) pointwise)

Version check_C28-1.0. Core oracle-1.0, constructions/strict_upper.py
(strict_upper-1.0). Every quantity is computed two ways and compared:

  (a) the closed-form formulas of the proof (Lemma R2 for S, Lemma R3 for F,
      the explicit shift rule of the Theorem for which c to use), evaluated
      as plain integer arithmetic;
  (b) `constructions.strict_upper.shift_word`, which builds and *executes*
      the actual word and asserts it sorts pi_h (ground truth, independent of
      the closed-form claims).

This is a proof-following checker (PLAN 9.11): it does not re-derive the
lemmas, it certifies that the closed forms match the executed construction
on every n, h tested, and that the resulting bound holds. A finite range of
n does not prove the "for all n" theorem; the algebraic argument in the
proof file is what does that. Coverage: exhaustive over h for each n in the
tested range of n (not sampled), both routes (H and carrier R).
"""

import os
import sys
import json
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "constructions"))
sys.path.insert(0, os.path.join(ROOT, "oracle"))

import strict_upper as su  # noqa: E402

CHECK_VERSION = "check_C28-1.0"
OUT_DIR = os.path.join(ROOT, "data", "runs", "check_C28")


def B(n):
    return n * (n - 1) // 2


def P(n):
    return n * n // 4


def refl(h, n):
    return tuple((h - i) % n for i in range(n))


def closed_form_S(hp, n):
    """Lemma R2."""
    if n % 2 == 1:
        return 3 * (n - 1) // 2
    return 3 * (n // 2 - 1) if hp % 2 == 0 else 3 * (n // 2)


def closed_form_F(hp, n):
    """Lemma R3."""
    if n % 2 == 1:
        return P(n)
    if n % 4 == 0:
        return P(n)
    # n = 2 mod 4
    return P(n) - 1 if hp % 2 == 0 else P(n) + 1


def H_of(c, n):
    c %= n
    return n if c == 0 else n + min(c, n - c) - 2


def chosen_shift(h, n):
    """Theorem's explicit shift rule: returns c (not h') to use for pi_h."""
    n_mod4 = n % 4
    if n % 2 == 1:
        return 1
    good_parity_is_odd = (n_mod4 == 0)   # n = 0 mod 4: odd h' is better; n = 2 mod 4: even h' is better
    h_even = (h % 2 == 0)
    if good_parity_is_odd:
        return 1 if h_even else 0
    else:
        return 0 if h_even else 1


def predicted_total(h, n):
    c = chosen_shift(h, n)
    hp = (h + c) % n
    F = closed_form_F(hp, n)
    S = closed_form_S(hp, n)
    return 2 * F - S + H_of(c, n), c, F, S


def run(n_max, log):
    out = {}
    for n in range(4, n_max + 1):
        t0 = time.time()
        rows = []
        max_H_excess = None
        max_R_excess = None
        for h in range(n):
            pi = refl(h, n)
            pred_val, c_pred, pred_F, pred_S = predicted_total(h, n)
            # ground truth at the predicted shift, via the executed construction
            sw_H = su.shift_word(pi, c_pred, "n1")
            st_H = sw_H["stats"]
            assert st_H["F_c"] == pred_F, (n, h, "F mismatch", st_H["F_c"], pred_F)
            assert st_H["S_c"] == pred_S, (n, h, "S mismatch", st_H["S_c"], pred_S)
            assert st_H["H"] == H_of(c_pred, n), (n, h, "H mismatch")
            actual_val = 2 * st_H["F_c"] - st_H["S_c"] + st_H["H"]
            assert actual_val == pred_val, (n, h, "total mismatch", actual_val, pred_val)
            assert actual_val <= B(n) + 1, (n, h, "BOUND FAILS (H-route)", actual_val, B(n))
            # true minimum over all c of the H-route quantity (must be <= our predicted value)
            true_min_H = min(2 * su.shift_word(pi, c, "n1")["stats"]["F_c"]
                              - su.shift_word(pi, c, "n1")["stats"]["S_c"]
                              + H_of(c, n) for c in range(n))
            assert true_min_H <= pred_val, (n, h, "predicted shift is not even a valid witness")
            h_exc = true_min_H - B(n)
            max_H_excess = h_exc if max_H_excess is None else max(max_H_excess, h_exc)
            # carrier route: R_c <= H(c) pointwise (trivial), so use it directly at c_pred too,
            # and also take the true min over all c for the carrier route (ground truth, executed)
            sw_R = su.shift_word(pi, c_pred, "carrier")
            R_at_pred = sw_R["stats"]["route_len"]
            assert R_at_pred <= st_H["H"], (n, h, "R_c > H(c) -- violates trivial Lemma R5", R_at_pred, st_H["H"])
            val_R_at_pred = 2 * st_H["F_c"] - st_H["S_c"] + R_at_pred
            assert val_R_at_pred <= B(n) + 1, (n, h, "BOUND FAILS (carrier route at predicted c)", val_R_at_pred, B(n))
            true_min_R = min(2 * su.shift_word(pi, c, "carrier")["stats"]["F_c"]
                              - su.shift_word(pi, c, "carrier")["stats"]["S_c"]
                              + su.shift_word(pi, c, "carrier")["stats"]["route_len"] for c in range(n))
            assert true_min_R <= B(n) + 1, (n, h, "BOUND FAILS (true min carrier route)", true_min_R, B(n))
            r_exc = true_min_R - B(n)
            max_R_excess = r_exc if max_R_excess is None else max(max_R_excess, r_exc)
            rows.append({"h": h, "c_pred": c_pred, "H_excess_true_min": h_exc, "R_excess_true_min": r_exc})
        out[n] = {"n": n, "B_n": B(n), "max_H_excess_true_min": max_H_excess,
                  "max_R_excess_true_min": max_R_excess, "rows": rows,
                  "seconds": round(time.time() - t0, 2)}
        log(f"n={n}: all {n} reflections OK; max true-min excess H-route {max_H_excess}, "
            f"carrier-route {max_R_excess}; predicted-shift formulas matched exactly ({out[n]['seconds']} s)")
    return out


def main():
    n_max = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    os.makedirs(OUT_DIR, exist_ok=True)
    logs = []

    def log(msg):
        print(msg, flush=True)
        logs.append(msg)

    t0 = time.time()
    res = run(n_max, log)
    meta = {"check": CHECK_VERSION, "construction": su.CONSTRUCTION_VERSION,
            "command": " ".join(sys.argv), "n_max": n_max,
            "seconds_total": round(time.time() - t0, 2)}
    with open(os.path.join(OUT_DIR, "report.json"), "w") as f:
        json.dump({"meta": meta, "results": res, "log": logs}, f, indent=1)
    log(f"PASS: closed-form formulas (Lemmas R2, R3) and the theorem's explicit shift rule "
        f"match the executed construction exactly for all reflections, 4<=n<={n_max}; "
        f"bound <= B_n+1 holds for both routes at every n, h tested. Total {round(time.time()-t0,2)} s")


if __name__ == "__main__":
    main()
