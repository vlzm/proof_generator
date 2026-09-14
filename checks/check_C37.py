"""check_C37.py -- independent checker for C37 (docs/proofs/C37_inversion_spearman.md).

Everything is recomputed from the definitions: inversions of the line are counted
by the naive double loop (no incremental formulas), V and D are summed directly.
All arithmetic is integer: V4 = 4*V, Q4 = 4*Q.

Parts
  A  identity (C37.1)     2n*inv(a,b) + V4(a,b) is independent of (a,b)
  B  corollaries (C37.3), (C37.4), the lag decomposition of section 6 and the
     discrepancy form (C37.3b) of section 6a
  C  constant (C37.5), (C37.6)
  D  equivalence (C37.7) and sufficiency (C37.8) against brute-force I(pi)
  E  reflections: closed forms of section 9 and equality in (C37.7), 4 <= n <= 60
  F  H13-I itself by brute force, 4 <= n <= 8, and cross-check of
     experiments/line_imax.c (compiled and run here) for 4 <= n <= 10

Usage: python3 checks/check_C37.py [--nmax 8] [--reflmax 60] [--skip-c]
Version check_C37-1.0.
"""

import argparse
import itertools
import json
import os
import subprocess
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "data", "runs", "check_C37")
VERSION = "check_C37-1.0"


# ---------- primitives, all independent of experiments/ ----------

def line(pi, a, b):
    """s_i and t_i for every point i (i is the index of the loop)."""
    n = len(pi)
    s = [(i - a) % n for i in range(n)]
    t = [(pi[i] - b) % n for i in range(n)]
    return s, t


def inv_naive(pi, a, b):
    n = len(pi)
    s, t = line(pi, a, b)
    c = 0
    for i in range(n):
        for j in range(n):
            if s[i] < s[j] and t[i] > t[j]:
                c += 1
    return c


def V4(pi, a, b):
    """4 * V(a,b) = sum_i (2 s_i - (n-1)) (2 t_i - (n-1))."""
    n = len(pi)
    s, t = line(pi, a, b)
    return sum((2 * s[i] - (n - 1)) * (2 * t[i] - (n - 1)) for i in range(n))


def Dsq(pi, a, b):
    n = len(pi)
    s, t = line(pi, a, b)
    return sum((s[i] - t[i]) ** 2 for i in range(n))


def reps(n):
    """One representative per (position rotation, value shift) orbit is enough for
    statements invariant under the torus action; we use all pi with pi(0) = 0,
    which is a superset of a transversal."""
    for tail in itertools.permutations(range(1, n)):
        yield (0,) + tail


def floor_bound(n):
    return ((n - 1) ** 2) // 4


def c_n_times_4n2(n):
    """4 n^2 c_n with c_n = (n^2-1)/4 - floor((n-1)^2/4)  (integer)."""
    return n * n * (n * n - 1) - 4 * n * n * floor_bound(n)


# ---------- parts ----------

def part_ABCD(nmax, log):
    stats = []
    for n in range(4, nmax + 1):
        t0 = time.time()
        cnt = 0
        for pi in reps(n):
            cnt += 1
            # --- A: identity
            ref = None
            Vmax4 = None
            invs = []
            Vs = []
            Ds = []
            for a in range(n):
                for b in range(n):
                    iv = inv_naive(pi, a, b)
                    v4 = V4(pi, a, b)
                    d = Dsq(pi, a, b)
                    invs.append(iv)
                    Vs.append(v4)
                    Ds.append(d)
                    val = 2 * n * iv + v4
                    if ref is None:
                        ref = val
                    elif val != ref:
                        raise AssertionError(f"A failed n={n} pi={pi} cut=({a},{b})")
                    if Vmax4 is None or v4 > Vmax4:
                        Vmax4 = v4
                    # --- B: (C37.3) 2D = n(n^2-1)/3 - V4
                    assert 2 * d == n * (n * n - 1) // 3 - v4, f"B1 n={n} pi={pi}"
            # --- B: (C37.4) n*inv - D constant
            nconst = None
            for k in range(n * n):
                val = n * invs[k] - Ds[k]
                if nconst is None:
                    nconst = val
                elif val != nconst:
                    raise AssertionError(f"B2 failed n={n} pi={pi}")
            # --- B3: lag decomposition of section 6 of the proof
            for a in range(n):
                for b in range(n):
                    c = (a - b) % n
                    delta = [((i - pi[i]) - c) % n for i in range(n)]
                    S = [i for i in range(n) if (i - a) % n < delta[i]]
                    assert len(S) * n == sum(delta), f"B3 |S| n={n} pi={pi}"
                    got = sum(d * d for d in delta) + n * sum(n - 2 * delta[i] for i in S)
                    assert got == Dsq(pi, a, b), f"B3 D n={n} pi={pi} cut=({a},{b})"
            # --- B4: discrepancy form (C37.3b), in 4n-scaled integers
            Nbox = [[sum(1 for x in range(A) if pi[x] < B) for B in range(n + 1)]
                    for A in range(n + 1)]
            # n*D(a,b) = n*N - a*b  (integer)
            nD = lambda A, B: n * Nbox[A][B] - A * B
            T = sum(i * pi[i] for i in range(n))
            gamma4 = 4 * T - 2 * n * (n - 1) * (n - 1) + n * (n - 1) ** 2   # 4*gamma
            for a in range(n):
                for b in range(n):
                    # 4*n*Phi = 4*(n*nD(a,b) - sum_b' nD(a,b') - sum_a' nD(a',b))
                    phi4n = 4 * (n * nD(a, b)
                                 - sum(nD(a, bb) for bb in range(1, n))
                                 - sum(nD(aa, b) for aa in range(1, n)))
                    assert V4(pi, a, b) == phi4n + gamma4, f"B4 n={n} pi={pi} cut=({a},{b})"
            # --- C: (C37.6)  12 n (n*const) = n^2(n^2-1) - 3 Q4
            Q4 = sum(V4(pi, i, pi[i]) for i in range(n))
            assert 12 * n * nconst == n * n * (n * n - 1) - 3 * Q4, f"C failed n={n} pi={pi}"
            # (C37.5): n^2 * inv_avg = n^2 (n^2-1)/4 - Q  ->  4 n^2 sum(inv) = n^2(n^2-1)n^2 ... use integers
            S = sum(invs)          # = n^2 * inv_avg
            assert 4 * S == n * n * (n * n - 1) - Q4, f"C2 failed n={n} pi={pi}"
            # --- D: (C37.7) equivalence and (C37.8) sufficiency
            I = min(invs)
            lhs = 2 * n * Vmax4 + Q4          # = 4n(2V* + Q/n)
            rhs = c_n_times_4n2(n)            # = 4n * (n c_n)
            assert (lhs >= rhs) == (I <= floor_bound(n)), f"D failed n={n} pi={pi}"
            # sufficiency (C37.8): 3 V* >= n c_n  =>  H13-I.
            # 3*V* >= n*c_n  <=>  3*V4* >= 4 n c_n  <=>  3 n V4* >= 4 n^2 c_n = rhs
            suff = 3 * n * Vmax4 >= rhs
            if suff:
                assert I <= floor_bound(n), f"D-suff failed n={n} pi={pi}"
        dt = time.time() - t0
        stats.append({"n": n, "representatives": cnt, "seconds": round(dt, 2)})
        log(f"  A-D ok for n={n}: {cnt} representatives, {dt:.1f}s")
    return stats


def part_E(reflmax, log):
    """Reflections: V(a,b) depends only on m = a+b, closed forms of section 9,
    and equality in (C37.7)."""
    rows = []
    for n in range(4, reflmax + 1):
        pi = [(-i) % n for i in range(n)]
        # direct check that V depends only on a+b, and the closed form
        # V(m) = -C(m) - n psi(-m)   (in 4x integers: V4(m) = -4C(m) - 4 n psi(-m))
        psi4 = lambda t: 2 * ((t % n)) - (n - 1)          # = 2*psi(t)
        C4 = lambda m: sum(psi4(t) * psi4((t + m) % n) for t in range(n))  # = 4*C(m)
        for a in range(min(n, 7)):
            for b in range(min(n, 7)):
                m = (a + b) % n
                got = V4(pi, a, b)
                want = -C4(m) - 2 * n * psi4((-m) % n)
                assert got == want, f"E1 n={n} a={a} b={b}: {got} != {want}"
        V4s = [V4(pi, 0, m) for m in range(n)]
        Vstar4 = max(V4s)
        Q4 = sum(V4(pi, i, pi[i]) for i in range(n))
        assert Q4 == n * V4s[0], f"E2 n={n}"
        # closed forms: 4*V(0) = n(n-1)(5-n)/3 ; V*-V(0) = n(n-2)^2/8 (even) or n(n-1)(n-3)/8 (odd)
        assert 3 * V4s[0] == n * (n - 1) * (5 - n), f"E3 n={n}"
        gap4 = Vstar4 - V4s[0]
        want_gap4 = n * (n - 2) ** 2 // 2 if n % 2 == 0 else n * (n - 1) * (n - 3) // 2
        assert gap4 == want_gap4, f"E4 n={n}: {gap4} != {want_gap4}"
        # equality in (C37.7)
        lhs = 2 * n * Vstar4 + Q4
        rhs = c_n_times_4n2(n)
        assert lhs == rhs, f"E5 n={n}: {lhs} != {rhs}"
        rows.append(n)
    log(f"  E ok for reflections, 4 <= n <= {reflmax} (equality in C37.7 at every n)")
    return rows


def part_F(nmax_py, skip_c, log):
    """H13-I by brute force in Python (4..nmax_py) and via experiments/line_imax.c."""
    py = {}
    for n in range(4, nmax_py + 1):
        best = -1
        arg = None
        for pi in reps(n):
            I = min(inv_naive(pi, a, b) for a in range(n) for b in range(n))
            if I > best:
                best = I
                arg = pi
        py[n] = {"max_I": best, "bound": floor_bound(n), "argmax": list(arg)}
        assert best <= floor_bound(n), f"H13-I violated at n={n}"
        log(f"  F python: n={n} max I = {best} <= {floor_bound(n)}, argmax {arg}")
    cres = {}
    if not skip_c:
        src = os.path.join(ROOT, "experiments", "line_imax.c")
        exe = os.path.join(OUT, "line_imax")
        os.makedirs(OUT, exist_ok=True)
        subprocess.run(["gcc", "-O2", "-o", exe, src], check=True)
        for n in range(4, 11):
            out = subprocess.run([exe, str(n)], capture_output=True, text=True, check=True)
            obj = json.loads(out.stdout.strip().splitlines()[-1])
            cres[n] = obj
            assert obj["violations"] == 0, f"line_imax reports a violation at n={n}"
            if n in py:
                assert obj["max_I"] == py[n]["max_I"], f"C/python disagree at n={n}"
            log(f"  F line_imax: n={n} max I = {obj['max_I']} (bound {obj['bound_floor_(n-1)^2/4']}), "
                f"unique class: {obj['argmax_count'] == 1}, {obj['seconds']}s")
    return py, cres


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=8, help="exhaustive n for parts A-D and F(python)")
    ap.add_argument("--reflmax", type=int, default=60)
    ap.add_argument("--skip-c", action="store_true")
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    lines = []

    def log(s):
        print(s)
        lines.append(s)

    t0 = time.time()
    log(f"{VERSION}: checking C37 (docs/proofs/C37_inversion_spearman.md)")
    log("A-D: identity, corollaries, constant, equivalence with H13-I")
    stats = part_ABCD(args.nmax, log)
    log("E: reflections")
    part_E(args.reflmax, log)
    log("F: H13-I by brute force and via experiments/line_imax.c")
    py, cres = part_F(min(args.nmax, 8), args.skip_c, log)
    dt = time.time() - t0
    log(f"PASS in {dt:.1f}s")

    report = {
        "version": VERSION,
        "claim": "C37",
        "proof": "docs/proofs/C37_inversion_spearman.md",
        "coverage": {
            "parts_A_D": f"all pi with pi(0)=0, all n^2 cuts, 4 <= n <= {args.nmax}",
            "part_E": f"reflections, 4 <= n <= {args.reflmax}",
            "part_F_python": f"4 <= n <= {min(args.nmax, 8)} exhaustive",
            "part_F_line_imax": "4 <= n <= 10 exhaustive" if not args.skip_c else "skipped",
        },
        "per_n": stats,
        "h13i_python": py,
        "h13i_line_imax": cres,
        "seconds": round(dt, 1),
        "result": "PASS",
    }
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(report, f, indent=1, ensure_ascii=False)
    with open(os.path.join(OUT, "log.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
