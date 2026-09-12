"""Finite checks of docs/proofs/C28_reflections.md (theorem C28: the carrier-route
variant of construction N1 sorts every reflection pi_h(i) = h - i with a word of
length <= B_n + 1, and <= B_n after free reduction; hence d(pi_h) <= B_n).

Every assert names the lemma / case of the proof it tests. The words are built
by constructions/strict_upper.py (strict_upper-1.0) and executed letter by
letter with the reference moves of oracle/moves.py (oracle-1.0). Passing
finite checks test the implementation and the formulas on the range; the
general argument is the text of the proof (AGENTS.md rule 5).

Usage: python3 checks/check_C28.py [--max-n 64] [--sigma-max-n 200] [--min-max-n 40]
Output: data/runs/check_C28/report.json and summary.md
"""

import argparse
import json
import os
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "oracle"))
sys.path.insert(0, os.path.join(ROOT, "constructions"))

from moves import apply_word, freely_reduce, identity, delta, sigma, CORE_VERSION  # noqa: E402
import strict_upper as su  # noqa: E402

CHECK_VERSION = "check_C28-1.0"
OUT_DIR = os.path.join(ROOT, "data", "runs", "check_C28")


def B(n):
    return n * (n - 1) // 2


def P(n):
    return n * n // 4


def reflection(h, n):
    return tuple((h - i) % n for i in range(n))


def chosen_shift(n, h):
    """Cases 1-6 of the proof: (case, c)."""
    m = n // 2
    if n % 2 == 1:
        return 1, 1
    if n % 4 == 0:
        return (2, 1) if h % 2 == 0 else (3, 0)
    if h % 2 == 1:
        return 4, 1
    if h % n == (m + 1) % n:
        return 6, 2
    return 5, 1


def predicted(n, hp):
    """Lemma 3 and Lemma 5: (fix, K, F_c, 2F_c - S_c) for f_c = pi_{h'}."""
    m = n // 2
    if n % 2 == 1:
        return 1, (n - 1) // 2, P(n), B(n) - (n - 1)
    if n % 4 == 0:
        if hp % 2 == 1:
            return 0, m, P(n), B(n) - n
        return 2, m - 1, P(n), B(n) - n + 3
    if hp % 2 == 1:
        return 0, m, P(n) + 1, B(n) - n + 2
    return 2, m - 1, P(n) - 1, B(n) - n + 1


def check_reflection(n, h, want_min):
    """Returns a dict of measured values; raises AssertionError on failure."""
    pi = reflection(h, n)
    m = n // 2
    case, c = chosen_shift(n, h)
    hp = (h + c) % n
    sw = su.shift_word(pi, c, "carrier")           # asserts: the word sorts pi
    st = sw["stats"]
    word = sw["word"]
    fix = sum(1 for i in range(n) if (2 * i - hp) % n == 0)
    K = st["n_cycles"]
    e_fix, e_K, e_F, e_2FS = predicted(n, hp)
    assert (fix, K) == (e_fix, e_K), ("Lemma 3 fix/K", n, h, fix, K, e_fix, e_K)
    assert st["F_c"] == sum(delta((2 * i) % n, hp, n) for i in range(n)), ("Lemma 5 identity F_c", n, h)
    assert st["F_c"] == e_F, ("Lemma 5 F_c", n, h, st["F_c"], e_F)
    assert st["S_c"] == 3 * K, ("Lemma 4/5 S_c = 3K", n, h)
    assert all(lw["cd"]["k"] == 2 and lw["cd"]["S"] == 3 for lw in sw["locals"]), ("Lemma 4: transpositions with S = 3", n, h)
    assert 2 * st["F_c"] - st["S_c"] == e_2FS, ("Lemma 5 table", n, h)
    assert st["local_len"] == e_2FS, ("(2.6) sum |W_C| = 2F - S", n, h)
    # Lemma 2: route bound
    assert st["route_len"] <= n - delta(0, c, n), ("Lemma 2 R_c <= n - delta(0,c)", n, h, c, st["route_len"])
    if case == 6:
        assert 1 not in st["carriers"], ("Case 6: vertex 1 is not a carrier", n, h, st["carriers"])
        assert st["theta_c"] == 0, ("Case 6: no antipodal pairs", n, h)
    bound = e_2FS + n - delta(0, c, n)
    expected_bound = {1: B(n), 2: B(n) - 1, 3: B(n), 4: B(n), 5: B(n) + 1, 6: B(n) - 1}[case]
    assert bound == expected_bound, ("Case arithmetic", n, h, case, bound, expected_bound)
    assert len(word) == st["len"] == 2 * st["F_c"] - st["S_c"] + st["route_len"] <= bound, ("length", n, h)
    red = freely_reduce(word)
    assert apply_word(pi, red) == identity(n), ("reduced word must sort", n, h)
    lr_junction = None
    if case == 5:
        anti = [lw for lw in sw["locals"] if lw["cd"]["antipodal_transposition"]]
        assert len(anti) == 1, ("Case 5: exactly one antipodal pair", n, h)
        a = anti[0]["carrier"]
        assert 0 <= a < m and a != 1, ("Case 5: antipodal carrier in [0,m), not 1", n, h, a)
        W = anti[0]["word"]
        assert W == "XL" * (m - 1) + "L" + "XL" * (m - 1), ("Lemma 6 word form", n, h, W)
        # the assembled word: W_a followed by a route letter R
        pos = word.find(W)
        assert pos >= 0 and pos + len(W) < len(word) and word[pos + len(W)] == "R", ("Case 5: W_a followed by R", n, h)
        lr_junction = pos + len(W) - 1
        assert len(red) <= len(word) - 2, ("Case 5: free reduction removes the LR pair", n, h)
    assert len(red) <= B(n), ("Theorem part 2: reduced length <= B_n", n, h, len(red))
    assert len(word) <= B(n) + 1, ("Theorem part 1: length <= B_n + 1", n, h, len(word))
    out = {"n": n, "h": h, "case": case, "c": c, "h'": hp, "fix": fix, "K": K, "F_c": st["F_c"], "S_c": st["S_c"],
           "R_c": st["route_len"], "len": len(word), "red": len(red), "len_minus_B": len(word) - B(n),
           "red_minus_B": len(red) - B(n), "N_X": st["N_X"], "N_rot": st["N_rot"], "LR_junction_index": lr_junction}
    if want_min:
        lens = []
        reds = []
        for cc in range(n):
            w = su.shift_word(pi, cc, "carrier")["word"]
            lens.append(len(w))
            reds.append(len(freely_reduce(w)))
        out["min_c_len_R_minus_B"] = min(lens) - B(n)
        out["min_c_red_R_minus_B"] = min(reds) - B(n)
        out["chosen_is_min_len"] = len(word) == min(lens)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=64)
    ap.add_argument("--sigma-max-n", type=int, default=200)
    ap.add_argument("--min-max-n", type=int, default=40)
    args = ap.parse_args()
    t0 = time.time()
    rows = []
    per_case = {}
    for n in range(4, args.max_n + 1):
        for h in range(n):
            r = check_reflection(n, h, want_min=(n <= args.min_max_n))
            rows.append(r)
            per_case.setdefault(r["case"], []).append(r["len_minus_B"])
    # sigma_n up to sigma_max_n: length <= B_n, N_X = floor((n-1)^2/4) (Lemma 8)
    sig = []
    for n in range(4, args.sigma_max_n + 1):
        r = check_reflection(n, 1, want_min=False)
        assert reflection(1, n) == sigma(n)
        assert r["len"] <= B(n) and r["red"] <= B(n), ("sigma_n: length <= B_n", n)
        assert r["N_X"] == (n - 1) ** 2 // 4, ("Lemma 8: N_X of sigma_n word", n, r["N_X"])
        assert r["N_rot"] <= n * n // 4, ("sigma_n: N_rot <= floor(n^2/4)", n, r["N_rot"])
        sig.append({"n": n, "len_minus_B": r["len_minus_B"], "red_minus_B": r["red_minus_B"], "N_X": r["N_X"], "N_rot": r["N_rot"], "c": r["c"]})
    elapsed = time.time() - t0
    summary = {
        "check": CHECK_VERSION, "core": CORE_VERSION, "construction": su.CONSTRUCTION_VERSION,
        "command": " ".join(sys.argv), "elapsed_s": round(elapsed, 1),
        "coverage": {
            "reflections_all_h": f"4<=n<={args.max_n}, all h ({len(rows)} inputs), chosen shift, words built, executed and reduced",
            "min_over_all_c": f"4<=n<={args.min_max_n}, all h, all c (report only)",
            "sigma_n": f"4<=n<={args.sigma_max_n}",
        },
        "per_case_len_minus_B": {str(k): sorted(set(v)) for k, v in sorted(per_case.items())},
        "max_len_minus_B": max(r["len_minus_B"] for r in rows),
        "max_red_minus_B": max(r["red_minus_B"] for r in rows),
        "chosen_c_not_minimal_count": sum(1 for r in rows if "chosen_is_min_len" in r and not r["chosen_is_min_len"]),
        "max_min_c_len_R_minus_B": max(r["min_c_len_R_minus_B"] for r in rows if "min_c_len_R_minus_B" in r),
        "max_min_c_red_R_minus_B": max(r["min_c_red_R_minus_B"] for r in rows if "min_c_red_R_minus_B" in r),
        "sigma_len_minus_B_values": sorted({s["len_minus_B"] for s in sig}),
        "sigma_red_minus_B_values": sorted({s["red_minus_B"] for s in sig}),
        "result": "PASS",
    }
    os.makedirs(OUT_DIR, exist_ok=True)
    json.dump({"summary": summary, "reflections": rows, "sigma": sig}, open(os.path.join(OUT_DIR, "report.json"), "w"), indent=1)
    lines = [f"# check_C28 — {summary['result']}", "",
             f"Версии: {CHECK_VERSION}, {su.CONSTRUCTION_VERSION}, {CORE_VERSION}; команда `{summary['command']}`; {elapsed:.1f} с.", "",
             "Покрытие: " + "; ".join(f"{k}: {v}" for k, v in summary["coverage"].items()), "",
             "Все assert'ы лемм 2–6, 8 и случаев 1–6 прошли. Значения `len − B_n` по случаям (выбранный сдвиг):", ""]
    for k, v in summary["per_case_len_minus_B"].items():
        lines.append(f"- случай {k}: {v}")
    lines += ["", f"Максимум `len − B_n` = {summary['max_len_minus_B']}, максимум `red − B_n` = {summary['max_red_minus_B']}.",
              f"Минимум по всем c при n <= {args.min_max_n}: max `min_c len_R − B_n` = {summary['max_min_c_len_R_minus_B']}, "
              f"max `min_c red_R − B_n` = {summary['max_min_c_red_R_minus_B']}; выбранный сдвиг не минимален по длине на "
              f"{summary['chosen_c_not_minimal_count']} входах (это не влияет на теорему).",
              f"sigma_n при 4 <= n <= {args.sigma_max_n}: `len − B_n` принимает значения {summary['sigma_len_minus_B_values']}, "
              f"N_X = floor((n−1)^2/4) на всех n."]
    open(os.path.join(OUT_DIR, "summary.md"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
