"""Independent finite checks of N1 (docs/incoming/LRX_STRICT_PROOF_RU-1.md)
for claims C16-C24 of CLAIMS.md (PLAN.md §9.10a).

This is a NEW implementation written from the proof text. It is not the
author's ``verify_lrx.py`` (not supplied with the incoming files). Every
assert is tagged with the N1 statement it checks. Passing all asserts is a
finite check of the implementation and a bug hunt; the general proof is
audited in docs/notes/strict_proof_audit.md.

Modes (can be combined):
  --cycles N        all nontrivial cycles on Z_m, 4<=m<=N, up to rotation of
                    the cycle notation: N1 (1.1), F<=nM, support connectivity,
                    (2.1)-(2.6), (2.5) action, (7.1), (7.2) classification
  --perms N         all permutations 4<=m<=N, all shifts c: words sort,
                    (3.3), averages (3.4), (4.4), (4.5), (5.1) at the
                    proof's h, (6.1) for all h, (6.4), (7.3), (8.1)-(8.3),
                    (T) with exact arithmetic; comparison with certified
                    tables (data/tables) when available (loss statistics)
  --sampled         random / perturbed-reflection permutations at larger n
                    (words built and executed for n<=31; only statistics and
                    lemmas for n=50,100) -- SAMPLED, not VERIFIED
  --scalar N        exact integer rounding checks of (T), (T') and §9 for
                    4<=n<=N (no floats)
  --reflections N   C24: s = 3(n-1)/2 on all reflections, 4<=n<=N
  --c22 N           counterexample search for C22 (2mu + E_bar >= n-1)
  --c23 B           C23: family of N1 §10.1 for n = 8b, 1<=b<=B
  --report PATH     write a JSON report with coverage fields
"""

import argparse
import itertools
import json
import os
import random
import sys
import time
from fractions import Fraction
from math import isqrt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "oracle"))
sys.path.insert(0, os.path.join(HERE, "..", "constructions"))
sys.path.insert(0, os.path.join(HERE, "..", "exact"))
sys.path.insert(0, os.path.join(HERE, "..", "bounds"))
from moves import apply_word, delta, CORE_VERSION  # noqa: E402
import strict_upper as SU  # noqa: E402
from known import strict_upper_bound, strict_upper_bound_simple, \
    target_diameter, ceil_sqrt  # noqa: E402


# --------------------------------------------------------------- helpers ---

def rhs_71(st, n):
    """Right-hand side of N1 (7.1) as an exact Fraction."""
    return (Fraction(4, 3) * st["M"] + Fraction(st["F"], 3 * n)
            + Fraction(7, 9) * st["E"]
            - (Fraction(2, 9) if st["antipodal_transposition"] else 0))


def sqrt_le(y, z):
    """Exact test  y <= sqrt(z)  for rational y and integer z >= 0."""
    if y <= 0:
        return True
    return y * y <= z


def sqrt_ge(y, z):
    """Exact test  y >= sqrt(z)  for rational y and integer z >= 0."""
    if y < 0:
        return False
    return y * y >= z


def r_n_radicand(n):
    P = n * n // 4
    return 3 * n * n - 4 * n + 1 - 4 * P


def s_lower_bound_83_holds(s, n):
    """Exact test of (8.3): s >= 5P/(3n) + 7 r_n / 9, where
    r_n = (sqrt(z) - (n-1))/2, z = 3n^2-4n+1-4P."""
    P = n * n // 4
    y = Fraction(9, 7) * (s - Fraction(5 * P, 3 * n))   # need y >= r_n
    return sqrt_ge(2 * y + (n - 1), r_n_radicand(n))


def A_of(x, n):
    P = n * n // 4
    return Fraction(5 * P, 3 * n) + Fraction(7, 9) * x


def B_of(x, n):
    P = n * n // 4
    return (Fraction(2 * (n - 1), 3) + Fraction(P, 3 * n)
            - (Fraction(5, 9) - Fraction(4, 3 * n)) * x
            - Fraction(4, 3 * n) * x * x)


def R_n_upper_int(n):
    """U_n = floor(R_n) via bounds/known.py (exact)."""
    return strict_upper_bound(n)


def load_table(n):
    tables = os.path.join(HERE, "..", "data", "tables")
    binp = os.path.join(tables, f"dist_n{n}.bin")
    metap = os.path.join(tables, f"dist_n{n}.json")
    if not (os.path.exists(binp) and os.path.exists(metap)):
        return None
    with open(metap) as f:
        meta = json.load(f)
    if not meta.get("certified"):
        return None
    with open(binp, "rb") as f:
        return f.read()


# ---------------------------------------------------------------- cycles ---

def all_cycles(n):
    """All nontrivial cycles on Z_n up to rotation of the notation: the
    smallest element first, then every ordering of the rest."""
    for k in range(2, n + 1):
        for subset in itertools.combinations(range(n), k):
            first, rest = subset[0], subset[1:]
            for order in itertools.permutations(rest):
                yield (first,) + order


def check_cycles(max_n, report):
    t0 = time.time()
    total = 0
    per_n = {}
    class_counts = {}
    for n in range(4, max_n + 1):
        cnt = 0
        for C in all_cycles(n):
            cnt += 1
            st = SU.cycle_stats(C, n)
            k, ds, loads = st["k"], st["ds"], st["loads"]
            # N1 §1: steps in (-n/2, n/2], antipodal positive
            for d in ds:
                assert -n < 2 * d <= n and d != 0, "N1 §1 step range"
            # N1 (1.1) load-jump rule, at every vertex of the circle
            on_cycle = {a: j for j, a in enumerate(C)}
            for v in range(n):
                jump = loads[v] - loads[(v - 1) % n]
                if v not in on_cycle:
                    assert jump == 0, "(1.1) off-cycle vertex"
                else:
                    j = on_cycle[v]
                    din, dout = ds[j - 1], ds[j]
                    exp = 2 if din < 0 < dout else (-2 if din > 0 > dout else 0)
                    assert jump == exp, "(1.1) on-cycle vertex"
            assert st["t"] == st["t_minus"], "N1 §1: #(-+) == #(+-)"
            assert st["E"] >= 0, "N1 §1: E >= 0"
            assert st["F"] == sum(loads) and st["F"] <= n * st["M"], \
                "N1 §1: F = sum loads <= nM"
            assert SU.positive_support_is_one_arc(loads), \
                "N1 §1/§7.1: positive support connected"
            # N1 §7.1 classification of small loads
            M, E = st["M"], st["E"]
            if M == 1:
                assert st["t"] == 0 and E == k, "§7.1: M=1 => one sign"
                if k == 2:
                    assert st["antipodal_transposition"], \
                        "§7.1: M=1,k=2 => antipodal transposition"
            if M <= 2 and st["t"] > 0:
                assert M == 2 and set(loads) == {0, 2}, "§7.1 loads in {0,2}"
                assert st["t"] == 1 and k == E + 2, "(7.2)"
                if E == 0:
                    assert k == 2 and st["F"] <= n, "§7.1: M=2,E=0"
                if E == 1:
                    assert k == 3 and st["F"] <= n and sum(ds) == 0, \
                        "§7.1: M=2,E=1 triangle"
            # N1 (7.1)
            assert Fraction(st["S"]) >= rhs_71(st, n), "(7.1)"
            key = (M if M <= 3 else "3+", min(E, 3), st["antipodal_transposition"])
            class_counts[str(key)] = class_counts.get(str(key), 0) + 1
            # N1 §2: carrier, contraction, word, (2.5), (2.6)
            lw = SU.local_word(C, n)
            assert SU.local_word_action_ok(lw, n), "(2.5) action of W_C"
            assert len(lw["word"]) == 2 * st["F"] - st["S"], "(2.6)"
            assert len(lw["reduced"]) <= len(lw["word"])
            # carrier edge carries load M (§2.1)
            assert loads[lw["carrier"]] == M, "§2.1 carrier edge has load M"
        per_n[n] = cnt
        total += cnt
        print(f"cycles n={n}: {cnt} cycles OK")
    print(f"cycles total {total} in {time.time()-t0:.1f}s")
    report["cycles"] = {"range": f"4<=n<={max_n}", "count": total,
                        "per_n": per_n, "coverage": "exhaustive up to "
                        "rotation of the cycle notation",
                        "class_counts": class_counts,
                        "seconds": round(time.time() - t0, 1)}


# ----------------------------------------------------- one permutation ---

def check_perm(pi, build_words=True, table=None, fact=None):
    """All N1 checks for one permutation. Returns the summary dict."""
    n = len(pi)
    P = n * n // 4
    shifts = SU.all_shifts(pi)
    ident = tuple(range(n))
    for s in shifts:
        if build_words:
            assert apply_word(pi, s["word"]) == ident, "(3.2)/(3.3): word sorts"
            assert apply_word(pi, s["reduced"]) == ident
        assert s["local_len"] == s["local_bound"], "(2.6) summed"
        assert s["route_len"] == s["H"], "§3 route length H(c)"
        assert s["len"] == s["bound"], "(3.3) attained by the raw word"
        assert s["reduced_len"] <= s["len"]
    sm = SU.perm_summary(pi, shifts)
    # (3.4) averages
    assert sm["F_bar"] == P, "§3: F_bar = P"
    assert sm["H_bar"] == n - 2 + Fraction(P + 2, n), "§3: H_bar"
    assert sm["avg_bound"] == 2 * P + n - 2 + Fraction(P + 2, n) - sm["s"], \
        "(3.4)"
    # §4: calE = V + 2B, B = sum theta_c, x = E_bar - 2 theta_bar >= 0
    assert sm["calE"] == sm["V"] + 2 * sm["B"], "(4.4)"
    assert sm["B"] == sum(s["theta_c"] for s in shifts), "(4.5) B = sum theta"
    assert sm["x"] == sm["E_bar"] - 2 * sm["theta_bar"] >= 0, "(4.5)"
    if n % 2:
        assert sm["B"] == 0 and sm["calE"] == sm["V"], "§4 odd n: calE = V"
    # §5: h = a + pi(a) at a = argmin e_i
    e = sm["e"]
    m = min(e)
    a = e.index(m)
    h = (a + pi[a]) % n
    V, x = sm["V"], sm["x"]
    T_h = sm["T"][h]
    if m == 0:
        assert 2 * T_h <= V, "(5.5)"
    else:
        assert T_h <= V + m * (m - 1), "(5.8)"
    assert T_h <= (1 - Fraction(1, n)) * V + Fraction(V * V, n * n), "(5.1)"
    # (5.2)/(5.3) pairwise, and (5.6) for exceptional partners of a
    tau = sm and [(-v) % n for v in pi]
    u = [(tau[i] - i) % n for i in range(n)]
    Q = V - n * m
    exceptional_partners = 0
    for i in range(n):
        for j in range(i + 1, n):
            assert 2 * abs(delta(i, j, n) - delta(tau[i], tau[j], n)) \
                <= e[i] + e[j], "(5.2)"
            si, sj = SU.signed_step(i, j, n), SU.signed_step(tau[i], tau[j], n)
            anti = (n % 2 == 0) and (abs(si) == n // 2 or abs(sj) == n // 2)
            exceptional = (si > 0) != (sj > 0) and not anti
            if not exceptional:
                assert delta(u[i], u[j], n) == abs(delta(i, j, n)
                                                   - delta(tau[i], tau[j], n)), \
                    "§5.1 non-exceptional equality"
                assert 2 * delta(u[i], u[j], n) <= e[i] + e[j], "(5.3)"
            elif a in (i, j) and m > 0:
                jj = j if a == i else i
                exceptional_partners += 1
                assert delta(u[jj], u[a], n) <= 2 * m + Fraction(e[jj] - m, 2) \
                    + Fraction(Q, 2 * m), "(5.6)"
    assert exceptional_partners <= m, "§5.4: b <= |B_a| = m"
    # §6: (6.1) for every h; (6.4)
    M_all = []
    for s in shifts:
        loads = [0] * n
        for i in range(n):
            j = (pi[i] + s["c"]) % n
            for ed in SU.arc_edges(i, SU.signed_step(i, j, n), n):
                loads[ed] += 1
        M_all.append(max(loads))
        assert s["sumM"] >= max(loads), "§6: sum of maxima >= max of sum"
    for hh in range(n):
        assert 2 * sum(M_all) >= n * (n - 1) - 2 * sm["T"][hh], "(6.1)"
    mu = sm["mu"]
    assert mu >= Fraction(P, n), "(6.4) first"
    assert mu >= Fraction(n - 1, 2) - Fraction(T_h, n), "(6.4) second"
    # (7.3), (8.1)-(8.3), (T)
    s_ = sm["s"]
    assert s_ >= Fraction(4, 3) * mu + Fraction(P, 3 * n) \
        + Fraction(7, 9) * x + Fraction(4, 3) * sm["theta_bar"], "(7.3)"
    assert s_ >= A_of(x, n), "(8.1)"
    assert s_ >= B_of(x, n), "(8.2)"
    assert s_lower_bound_83_holds(s_, n), "(8.3)"
    U = R_n_upper_int(n)
    assert sm["min_bound"] <= sm["avg_bound"], "min <= mean over c"
    assert sm["min_bound"] <= U, "(T): some shift word is within U_n"
    if build_words:
        assert sm["min_len"] <= U
    # C22 as a separate record, not an assert
    sm["c22_holds"] = (2 * mu + sm["E_bar"] >= n - 1)
    sm["h_proof"] = h
    sm["T_h_proof"] = T_h
    sm["min_error"] = m
    if table is not None:
        from bfs import rank_perm
        d = table[rank_perm(tuple(pi), fact)]
        assert d <= sm["min_reduced_len"] <= sm["min_len"], \
            "certified distance must not exceed any constructed word"
        sm["d_exact"] = d
    return sm


def check_perms(max_n, report, want_table=True):
    from bfs import factorials
    out = {}
    for n in range(4, max_n + 1):
        t0 = time.time()
        table = load_table(n) if want_table else None
        fact = factorials(n) if table else None
        U = R_n_upper_int(n)
        Bn = target_diameter(n)
        worst_len = -1
        worst_examples = []
        loss = {"len_minus_d": {}, "reduced_minus_d": {},
                "minbound_minus_avgbound": {}, "d_minus_minlen": {}}
        c22_fail = 0
        count = 0
        slack_83 = None
        for pi in itertools.permutations(range(n)):
            sm = check_perm(pi, True, table, fact)
            count += 1
            if not sm["c22_holds"]:
                c22_fail += 1
            ml = sm["min_len"]
            if ml > worst_len:
                worst_len, worst_examples = ml, [pi]
            elif ml == worst_len and len(worst_examples) < 5:
                worst_examples.append(pi)
            if table is not None:
                d = sm["d_exact"]
                for key, val in (("len_minus_d", ml - d),
                                 ("reduced_minus_d", sm["min_reduced_len"] - d)):
                    loss[key][val] = loss[key].get(val, 0) + 1
            gap = sm["avg_bound"] - sm["min_bound"]
            loss["minbound_minus_avgbound"][str(-gap)] = \
                loss["minbound_minus_avgbound"].get(str(-gap), 0) + 1
        rec = {
            "n": n, "permutations": count, "shifts_per_perm": n,
            "words_built": count * n, "U_n": U, "B_n": Bn,
            "max_min_c_len": worst_len, "worst_examples": worst_examples,
            "c22_violations": c22_fail,
            "loss_min_c_len_minus_d": {str(k): v for k, v in sorted(loss["len_minus_d"].items())},
            "loss_min_c_reduced_minus_d": {str(k): v for k, v in sorted(loss["reduced_minus_d"].items())},
            "table_used": table is not None,
            "seconds": round(time.time() - t0, 1),
        }
        out[n] = rec
        print(f"perms n={n}: {count} perms x {n} shifts OK; "
              f"max_pi min_c len = {worst_len} (B_n={Bn}, U_n={U}); "
              f"C22 violations={c22_fail}; table={table is not None}; "
              f"{rec['seconds']}s")
        if table is not None:
            print(f"   loss min_c len - d(pi): {rec['loss_min_c_len_minus_d']}")
            print(f"   loss min_c reduced - d(pi): "
                  f"{rec['loss_min_c_reduced_minus_d']}")
    report["perms"] = {"range": f"4<=n<={max_n}", "coverage": "exhaustive",
                       "per_n": out}


# --------------------------------------------------------------- sampled ---

def perturbed_reflection(n, rng, swaps):
    h = rng.randrange(n)
    pi = [(h - i) % n for i in range(n)]
    for _ in range(swaps):
        i, j = rng.randrange(n), rng.randrange(n)
        pi[i], pi[j] = pi[j], pi[i]
    return tuple(pi)


def check_sampled(report, seed=20260912, cases=20):
    rng = random.Random(seed)
    orders = [8, 9, 10, 15, 20, 31, 50, 100]
    out = {}
    for n in orders:
        t0 = time.time()
        build = n <= 31
        c22_fail = 0
        worst = -1
        for idx in range(cases):
            if idx % 2 == 0:
                pi = tuple(rng.sample(range(n), n))
            else:
                pi = perturbed_reflection(n, rng, rng.randrange(0, 4))
            sm = check_perm(pi, build_words=build)
            if not sm["c22_holds"]:
                c22_fail += 1
            worst = max(worst, sm["min_len"] - target_diameter(n))
        out[n] = {"cases": cases, "words_executed": build,
                  "max_excess_over_B_n_of_min_c_len": worst,
                  "c22_violations": c22_fail,
                  "seconds": round(time.time() - t0, 1)}
        print(f"sampled n={n}: {cases} cases OK (words executed={build}); "
              f"max (min_c len - B_n) = {worst}; C22 violations={c22_fail}; "
              f"{out[n]['seconds']}s")
    report["sampled"] = {"seed": seed, "coverage": "SAMPLED", "per_n": out}


# ---------------------------------------------------------------- scalar ---

def check_scalar(max_n, report):
    """Exact integer checks of (T), (T') and N1 §9 for 4 <= n <= max_n."""
    t0 = time.time()
    for n in range(4, max_n + 1):
        P = n * n // 4
        z = r_n_radicand(n)
        # radicand positive, r_n >= 0 (§8)
        assert z > 0 and z >= (n - 1) ** 2, "§8: real nonnegative r_n"
        # r_n solves r^2 + (n-1) r = n(n-1)/2 - P  (checked symbolically:
        # z = (n-1)^2 + 4(n(n-1)/2 - P))
        assert z == (n - 1) ** 2 + 4 * (n * (n - 1) // 2 - P), "§8 root"
        U = strict_upper_bound(n)
        Tp = strict_upper_bound_simple(n)
        # (T') from (T): U_n <= B_n + floor(kappa n) - (2 odd | 1 even)
        assert U <= Tp, "(T') must follow from (T)"
        # §9: for odd n >= 7 the slack is strictly more than 1 below
        # the next integer? No: only U <= Tp is claimed; we also record
        # the exact difference Tp - U.
        if n == 5:
            assert U == 13 and Tp == 13, "§9: n = 5"
        # sanity: T' vs known formulas at control values (PROBLEM §5.2)
    ctrl = {4: 9, 7: 27, 20: 211, 100: 5065}
    for n, v in ctrl.items():
        assert strict_upper_bound(n) == v
    diffs = {}
    for n in range(4, min(max_n, 60) + 1):
        diffs[n] = strict_upper_bound_simple(n) - strict_upper_bound(n)
    # §9 real inequalities, exactly: odd n>=7: R_n < B_n + kappa n - 2;
    # even n>=4: R_n < B_n + kappa n - 1.  With q = (sqrt2 - 1)/2:
    # odd:  R_n - B_n - kappa n = -5/2 + 7q/9 + 13/(6n)
    # even: R_n - B_n - kappa n = -2 + 7q/9 + 2/n + 7 eps_n/9
    # We verify the *integer* consequence directly above (U <= Tp) for all
    # n; here we additionally verify the two closed forms of §9 at a few n
    # by exact algebra: R_n*18n is a - 7n sqrt(z) (see bounds/known.py).
    for n in range(4, max_n + 1):
        P = n * n // 4
        a = 36 * n * P + 25 * n * n - 43 * n - 12 * P + 36
        z = r_n_radicand(n)
        if n % 2:
            assert z == 2 * (n - 1) ** 2, "§9 odd radicand"
        else:
            assert z == 2 * (n - 1) ** 2 - 1, "§9 even radicand"
        # R_n = (a - 7n sqrt z)/(18n); (T') integer:
        # U_n = floor(R_n) computed with ceil_sqrt (exact)
        assert strict_upper_bound(n) == (a - ceil_sqrt(49 * n * n * z)) // (18 * n)
    print(f"scalar 4<=n<={max_n}: (T), (T') exact rounding OK "
          f"in {time.time()-t0:.1f}s; T'-U for n<=60: {diffs}")
    report["scalar"] = {"range": f"4<=n<={max_n}", "coverage": "exhaustive",
                        "arithmetic": "integer isqrt, no floats",
                        "Tprime_minus_U_n_le_60": diffs}


# ----------------------------------------------------------- reflections ---

def check_reflections(max_n, report):
    t0 = time.time()
    for n in range(4, max_n + 1):
        P = n * n // 4
        Bn = target_diameter(n)
        for h in range(n):
            pi = tuple((h - i) % n for i in range(n))
            sm = SU.perm_summary(pi)
            assert sm["s"] == Fraction(3 * (n - 1), 2), "C24: s = 3(n-1)/2"
            expr = 2 * P + n - 2 + Fraction(P + 2, n) - sm["s"]
            if n % 2 == 0:
                assert expr == Bn + Fraction(P, n) - Fraction(1, 2) + Fraction(2, n), \
                    "C24 even closed form"
            else:
                assert expr == Bn + Fraction(P, n) - 1 + Fraction(2, n), \
                    "C24 odd closed form"
            # sigma_n (h = 1): also record min_c bound vs B_n
            if h == 1:
                sig_min = sm["min_bound"]
        print(f"reflections n={n}: all h OK; sigma_n min_c bound - B_n = "
              f"{sig_min - Bn}")
    report["reflections"] = {"range": f"4<=n<={max_n}", "coverage":
                             "all reflections pi(i)=h-i", "seconds":
                             round(time.time() - t0, 1)}


# ------------------------------------------------------------------- C22 ---

def check_c22(max_n, report, samples=200, seed=20260912):
    rng = random.Random(seed)
    out = {}
    for n in range(4, max_n + 1):
        fails = []
        tight = 0
        it = itertools.permutations(range(n)) if n <= 7 else \
            (tuple(rng.sample(range(n), n)) for _ in range(samples))
        cnt = 0
        for pi in it:
            cnt += 1
            sm = SU.perm_summary(pi)
            lhs = 2 * sm["mu"] + sm["E_bar"]
            if lhs < n - 1:
                fails.append((pi, str(lhs)))
            elif lhs == n - 1:
                tight += 1
        out[n] = {"checked": cnt, "coverage": "exhaustive" if n <= 7 else
                  f"SAMPLED {samples}", "violations": len(fails),
                  "equality_cases": tight, "first_violations": fails[:5]}
        print(f"C22 n={n}: {cnt} perms ({out[n]['coverage']}), "
              f"violations={len(fails)}, equality={tight}")
    report["c22"] = out


# ------------------------------------------------------------------- C23 ---

def check_c23(max_b, report):
    out = {}
    for b in range(1, max_b + 1):
        n = 8 * b
        tau = [0] * n
        for k in range(8):
            for j in range(b):
                tau[b * k + j] = b * ((3 * k) % 8) + j
        pi = tuple((-t) % n for t in tau)
        assert sorted(pi) == list(range(n))
        re = SU.reflection_errors(pi)
        assert re["tau"] == tau
        assert all(ei == 2 * b for ei in re["e"]), "C23: e_i = 2b"
        assert re["V"] == n * n // 4, "C23: V = n^2/4"
        Ts = [SU.T(pi, h) for h in range(n)]
        assert min(Ts) == n * n // 4, "C23: T_min = n^2/4"
        B = SU.antipodal_pairs_preserved(pi)
        assert B == n // 2, "C23: B = n/2"
        calE = re["V"] + 2 * B
        if n <= 24:
            shifts = SU.all_shifts(pi)
            assert sum(s["E_c"] for s in shifts) == calE, "C23: calE via cycles"
        assert calE == re["V"] + n, "C23: calE = V + n"
        assert Fraction(min(Ts), calE) == Fraction(n, n + 4), "C23 ratio"
        out[n] = {"V": re["V"], "T_min": min(Ts), "calE": calE,
                  "cycles_checked": n <= 24}
        print(f"C23 n={n}: V={re['V']}, T_min={min(Ts)}, calE={calE}, "
              f"ratio={Fraction(min(Ts), calE)}")
    report["c23"] = out


# ------------------------------------------------------------------ main ---

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cycles", type=int)
    ap.add_argument("--perms", type=int)
    ap.add_argument("--no-table", action="store_true")
    ap.add_argument("--sampled", action="store_true")
    ap.add_argument("--sampled-cases", type=int, default=20)
    ap.add_argument("--scalar", type=int)
    ap.add_argument("--reflections", type=int)
    ap.add_argument("--c22", type=int)
    ap.add_argument("--c23", type=int)
    ap.add_argument("--report")
    args = ap.parse_args()
    report = {"core_version": CORE_VERSION,
              "construction_version": SU.CONSTRUCTION_VERSION,
              "command": " ".join(sys.argv),
              "mandatory_example": apply_word((2, 0, 1), "XL") == (2, 1, 0)}
    assert report["mandatory_example"]
    if args.cycles:
        check_cycles(args.cycles, report)
    if args.perms:
        check_perms(args.perms, report, want_table=not args.no_table)
    if args.sampled:
        check_sampled(report, cases=args.sampled_cases)
    if args.scalar:
        check_scalar(args.scalar, report)
    if args.reflections:
        check_reflections(args.reflections, report)
    if args.c22:
        check_c22(args.c22, report)
    if args.c23:
        check_c23(args.c23, report)
    if args.report:
        os.makedirs(os.path.dirname(os.path.abspath(args.report)), exist_ok=True)
        with open(args.report, "w") as f:
            json.dump(report, f, indent=1, default=str)
        print(f"report written: {args.report}")
    print("PASS")


if __name__ == "__main__":
    main()
