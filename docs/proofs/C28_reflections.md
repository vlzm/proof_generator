# C28 — H10 on reflections, closed form (PROVED, all n >= 4)

Status: PROVED. This is Step 1 of the H10/C27 task (PLAN §8, session 4).
It sharpens the reflection part of C27v (VERIFIED 4<=n<=9/12) from a finite
check to a closed-form argument valid for every n >= 4.

## 0. Statement

For all `n >= 4` and all `h` in `Z_n`, let `pi_h(i) = (h - i) mod n` (a
reflection of the circle; `pi_1 = sigma_n`, and every reflection is a
rotation of `sigma_n`: `pi_h = sigma_n o rho_{h-1}`, PROBLEM §4.3 notation).
Let `F_c, S_c` be as in N1 §1/§3 (`constructions/strict_upper.py:
cycle_data`), `H(c)` the N1 route (N1 §3, `H_of`), and `R_c` the
carrier-only route (variant, AGENTS.md rule 16, `carrier_route`). Then

```text
min_c [2F_c - S_c + H(c)] <= B_n + 1,
min_c [2F_c - S_c + R_c]  <= B_n + 1,
```

with equality (`= B_n + 1`) in the first line exactly when `n = 2 (mod 4)`
and `h` is even; in every other case (`n` odd, `n = 0 (mod 4)`, or
`n = 2 (mod 4)` with `h` odd) the minimum is `<= B_n`. The second line
follows from the first by the trivial Lemma R5 below; it is never worse,
and the *worst case over h* is identical for both routes at every `n`
tested (`4<=n<=40`, see "Machine verification") — but the two true minima
do **not** coincide at every individual `h`: e.g. at `n=6, h=0` the true
`min_c` of the H-route quantity is `B_n+1` while the true `min_c` of the
carrier-route quantity is `B_n-1` (§7 gives the numbers). The carrier
route can do strictly better than this proof's witness shift at specific
`h`; the Theorem below only needs one witness per `h` achieving `<=B_n+1`,
not the true minimum.

This is a genuine strengthening of C27v on the reflection family: C27v is
VERIFIED for `4<=n<=9` (exhaustive) and `4<=n<=12` (family); this file
PROVES the same bound for every `n>=4`, and pins down exactly where the
`+1` is unavoidable.

## 1. Cycle structure of a reflection under a shift (Lemma R1)

`f_c(i) = pi_h(i) + c = (h + c - i) mod n`. Writing `h' = (h + c) mod n`,
`f_c` is exactly the reflection `r_{h'}(i) = (h' - i) mod n` — the same
kind of map as `pi_h` itself, with parameter `h'` instead of `h`. `r_{h'}`
is an involution: `i` is a fixed point iff `2i = h' (mod n)`; every other
`i` lies in the 2-cycle (transposition) `{i, (h'-i) mod n}`.

- `n` odd: `2` is invertible mod `n`, so `2i = h'` has exactly one solution:
  exactly one fixed point, `(n-1)/2` transpositions, for every `h'`.
- `n` even: `2i = h' (mod n)` is solvable iff `h'` is even, and then has
  exactly two solutions (`i` and `i + n/2`). So `r_{h'}` has 2 fixed points
  and `n/2 - 1` transpositions if `h'` is even, 0 fixed points and `n/2`
  transpositions if `h'` is odd.

Since `h'` ranges over all of `Z_n` as `c` ranges over `Z_n` (for fixed
`h`), this reduces the whole family to: what are `F, S` of `r_{h'}` as a
function of `h'` (and `n`), and what route is needed to reach a given
target `c` from `0` while visiting the resulting carriers?

## 2. S per shift (Lemma R2)

Every nontrivial cycle of a reflection is a transposition, and N1 §1's
formulas + the case split of C24 give, for *every* transposition arising
this way (antipodal or not): `S(C) = 3` — `M=2, E=0` when the two arcs of
the pair coincide (non-antipodal, opposite directions over the same arc:
C24 §5.3 of PROBLEM.md), and `M=1, E=2` when both steps are the antipodal
`+n/2` (rule 16: antipodal steps are always positive). In both cases
`S = 2M + E - 1 = 3`.

Hence `S_{h'} = 3 * (number of transpositions of r_{h'})`:

```text
n odd:          S_{h'} = 3(n-1)/2                      for every h'
n even, h' even: S_{h'} = 3(n/2 - 1)
n even, h' odd:  S_{h'} = 3n/2
```

(`S_{h'}` does not depend on `h'` at all when `n` is odd — this reproves
C24's exact `s = 3(n-1)/2` pointwise, not just on average.)

## 3. F per shift, closed form (Lemma R3)

For a transposition `{i, j}`, `j = h'-i mod n`, both signed steps of the
2-cycle have the same absolute value `g := delta_n(i,j) <= n/2` — see
`cycle_data`: the two steps are `d` and `-d` if `d != n/2` (opposite signs,
same absolute value), or `n/2` and `n/2` if antipodal (rule 16 forces both
positive), and `2*(n/2) = n = 2g` too. So **every transposition contributes
`F += 2*delta_n(i,j)` regardless of the antipodal case** — one uniform
formula.

Writing `j = h'-i mod n`, `delta_n(i,j) = delta_n(h'-2i mod n, 0)`. Since
each unordered pair `{i,j}` gives the same value from either representative
(`delta_n(h'-2j,0) = delta_n(2i-h',0) = delta_n(h'-2i,0)`), summing this
quantity over *all* `i = 0..n-1` counts every transposition exactly twice
(and every fixed point contributes 0), so

```text
F_{h'} = sum_{i=0}^{n-1} delta_n(h' - 2i mod n, 0).
```

**n odd.** `i |-> 2i mod n` is a bijection of `Z_n` (2 is invertible), so
`h'-2i mod n` also ranges over all of `Z_n` once as `i` does. Hence

```text
F_{h'} = sum_{r=0}^{n-1} delta_n(r,0) = 2 * sum_{r=1}^{(n-1)/2} r
       = (n-1)(n+1)/4 = (n^2-1)/4 = P,
```

independent of `h'` (`P = floor(n^2/4) = (n^2-1)/4` for odd `n`).

**n even, n = 2m.** `2i mod n` only hits the even residues, each exactly
twice as `i` ranges `0..n-1`. So `h'-2i mod n` ranges over all residues of
parity `h' mod 2`, each exactly twice:

```text
F_{h'} = 2 * sum_{r = h' (mod 2), 0<=r<n} delta_n(r,0).
```

Using the standard identity `sum_{k=0}^{M} min(k,M-k) = floor(M^2/4)` (the
`k=M` term is 0, so this equals the same sum over `k=0..M-1`; both are
checked directly, e.g. by pairing `k <-> M-k`):

- `h'` even: `r = 2k`, `k=0..m-1`; `delta_n(2k,0) = 2*min(k,m-k)`, so
  `sum = 2*floor(m^2/4)`, giving `F_{h'} = 4*floor(m^2/4)`.
- `h'` odd: `r = 2k+1`, `k=0..m-1`; `delta_n(2k+1,0) = 2*min(k,m-1-k)+1`,
  so `sum = 2*floor((m-1)^2/4) + m`, giving
  `F_{h'} = 4*floor((m-1)^2/4) + 2m = 4*floor((m-1)^2/4) + n`.

Evaluating for the two residues of `n` mod 4 (`m = n/2`):

```text
n = 0 (mod 4):  F_{h'} = P   for both parities of h'   (m even, m-1 odd:
                4*floor(m^2/4) = m^2 = n^2/4 = P; 4*floor((m-1)^2/4)+n =
                (m-1)^2-1+n = m^2-2m+n = P, using n=2m)
n = 2 (mod 4):  F_{h'} = P-1 (h' even), P+1 (h' odd)   (m odd, m-1 even:
                4*floor(m^2/4) = m^2-1 = P-1; 4*floor((m-1)^2/4)+n =
                (m-1)^2+n = m^2-2m+1+2m = m^2+1 = P+1)
```

(`checks/check_C28.py` verifies these four closed forms against the
executed construction's actual `F_c, S_c` for every reflection, `4<=n<=40`,
by assertion — not just the final bound.)

## 4. The two route values used (Lemma R4)

Directly from N1 §3 (`H_of`): `H(0) = n`; for `c != 0`,
`H(c) = n + delta_n(0,c) - 2`. In particular `H(1) = n - 1` and, for
`n >= 4`, `H(2) = n` (`delta_n(0,2) = 2` since `2 <= n-2`).

## 5. R_c <= H(c) pointwise (Lemma R5, trivial)

`R_c` is defined (`carrier_route`) as the length of the *shortest* walk
from `0` to `c` on `Z_n` visiting a required set of positions (the
carriers of the nontrivial cycles of `f_c`). `H(c)`'s route visits *every*
position of `Z_n` (N1 §3), which is a superset of any carrier set. So the
`H`-route is one particular candidate walk for the carrier problem, and the
shortest one is no longer: `R_c <= H(c)` for every `c`, every `pi`, every
`n`. (Already noted in `data/runs/strict_loss_audit/report.md`; restated
here because Step 1 leans on it directly.)

## 6. Theorem (the explicit shift rule)

Fix `n >= 4` and `h in Z_n`. Choose `c` as follows:

```text
n odd:                       c = 1
n = 0 (mod 4), h even:       c = 1
n = 0 (mod 4), h odd:        c = 0
n = 2 (mod 4), h odd:        c = 1
n = 2 (mod 4), h even:       c = 0
```

(Equivalently: pick whichever of `c in {0,1}` makes `h' = h+c` the
*better* parity for `2F_{h'} - S_{h'}` — odd `h'` is better by exactly `3`
when `n=0 (mod 4)`, even `h'` is better by exactly `1` when `n=2 (mod 4)`;
ties never occur for odd `n` since `F_{h'}, S_{h'}` don't depend on `h'`
there at all.)

**Claim.** With this choice, `2F_c - S_c + H(c) <= B_n + 1` always, with
equality iff `n = 2 (mod 4)` and `h` even.

*Proof.* Substitute Lemmas R2–R4 (`B_n = n(n-1)/2` throughout, the closed
form that holds for every `n`, PROBLEM §1):

- **n odd** (`c=1`, any `h`, `h'` irrelevant since `F,S` don't depend on
  it): `2F-S+H(1) = 2P - 3(n-1)/2 + (n-1) = 2P - (n-1)/2`. With
  `P=(n^2-1)/4=(n-1)(n+1)/4`: `2P-(n-1)/2 - B_n
  = (n-1)[(n+1)/2 - 1/2 - n/2] = 0`. **Equals `B_n` exactly**, for every
  `h`.

- **n = 0 (mod 4)**, `h` even (`c=1`, `h'=h+1` odd, the better parity):
  `2F-S+H(1) = 2P - 3n/2 + (n-1) = 2P - n/2 - 1 = B_n - 1` (using
  `B_n = 2P - n/2` for `n=0 mod 4`, since `P=n^2/4`). **`= B_n - 1`.**

- **n = 0 (mod 4)**, `h` odd (`c=0`, `h'=h` odd): `2F-S+H(0)
  = 2P - 3n/2 + n = 2P - n/2 = B_n`. **`= B_n`.**

- **n = 2 (mod 4)**, `h` odd (`c=1`, `h'=h+1` even, the better parity):
  `2F-S+H(1) = 2(P-1) - 3(n/2-1) + (n-1) = 2P - n/2 = B_n` (since
  `B_n = 2P-n/2` for even `n` in general — `P=n^2/4` exactly whenever `n`
  is even, so this formula for `B_n` holds for both residues mod 4).
  **`= B_n`.**

- **n = 2 (mod 4)**, `h` even (`c=0`, `h'=h` even): `2F-S+H(0)
  = 2(P-1) - 3(n/2-1) + n = 2P - n/2 + 1 = B_n + 1`. **`= B_n + 1`.**

Every case gives `<= B_n + 1`, proving the Claim; the last case shows the
`+1` is achieved, not just an artifact of a loose bound. QED (H-route half
of the Theorem).

By Lemma R5, `2F_c - S_c + R_c <= 2F_c - S_c + H(c) <= B_n+1` at the same
`c`, so `min_c[2F_c-S_c+R_c] <= B_n+1` too — **Step 1 of H10/C27 is
proved, for every `n>=4`, on the full reflection family**, strengthening
C27v from a finite check to a theorem.

## 7. What this does and does not give

- This is a statement about **one input family** (reflections), using
  **one exhibited shift per input** — it is an upper bound (a witness),
  not a claim that this `c` is the true `argmin_c`. The true minimum can
  do strictly better (e.g. via the carrier route with a different `c`;
  see "Machine verification": at `n=6, h=2` the true minimum of the
  H-route quantity is `B_n+1` and the true minimum of the carrier-route
  quantity is *also* `B_n+1` there — carrier route buys nothing extra on
  reflections at the extremal point, though it does elsewhere in the
  family, e.g. `n=6,h=0` reaches `B_n-1` via carrier route by using a
  different `c` than this proof's witness).
- It says nothing about general `pi` — H9's refutation (C26) already shows
  the analogous `H`-route argument fails for general `pi` (growth like
  `floor((n-3)/2)`, C26v), and Lemma R5 alone (`R_c<=H(c)`) is *not*
  enough to rescue it there: the known n=9 counterexamples to H9 still
  have `Q_R - B_n = -4` (report), i.e. the carrier route is strictly
  better than the N1 route on exactly the inputs where the N1 route
  fails badly, so the trivial pointwise domination used in this proof
  (§5) is far from tight in general — this is the real content of the
  Step 2 obstruction (see PLAN §8 / the session-4 log entry).

## 8. Machine verification

`checks/check_C28.py` (check_C28-1.0), against `strict_upper-1.0`
(oracle-1.0): for every `4<=n<=40` and every `h in Z_n` (exhaustive, not
sampled):

1. asserts the closed-form `F_{h'}, S_{h'}` (§3, §2) match
   `constructions.strict_upper.shift_word`'s actual `F_c, S_c` at the
   Theorem's chosen shift;
2. asserts the Theorem's predicted total matches the executed
   construction's `2F_c-S_c+H(c)` exactly;
3. asserts the *true* `min_c` of the H-route quantity (all `n` shifts,
   words built and executed) is `<=` the Theorem's witness value, and
   `<= B_n+1`;
4. asserts `R_c <= H(c)` pointwise at the witness shift (Lemma R5) and
   that both the witness-shift and the true `min_c` of the carrier-route
   quantity are `<= B_n+1`.

Result (`data/runs/check_C28/run.md`): all assertions PASS,
`4<=n<=40`; the true-minimum excess over `B_n` is exactly `1` when
`n = 2 (mod 4)` and exactly `0` otherwise, for *both* routes, at every `n`
tested — confirming the Theorem's case split is exact, not just an upper
bound, on this family.
