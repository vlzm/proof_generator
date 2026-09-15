# C37 — exact formula for the value-shift half of the double cut (session 9)

Date: 15.09.2026 (session 9). Core oracle-1.0 (unused directly; this is pure
combinatorics on permutations of `Z_n`, no board moves involved). Companion
to `docs/notes/h13_line_model.md` (H13-I attempt) and `PLAN.md` H13 row.

## 0. Statement

Notation as in `docs/notes/h13_line_model.md` §0 and PLAN's H13' row: `pi`
a permutation of `Z_n`; for a position cut `a in Z_n` and a value cut
`b in Z_n`,

```
w^{a,b}_j = (pi((a+j) mod n) - b) mod n,   j = 0..n-1
inv(a,b)  = inv_count(w^{a,b})
I(pi)     = min_{a,b} inv(a,b)
```

(`a = q+1`, `b = q+1-c` relates this to the `(q,c)` notation of the line
model.)

Fix `a`. Let `pinv = pi^{-1}` and `sigma_a(t) = (pinv(t) - a) mod n` for
`t = 0..n-1` (so `sigma_a` is the inverse permutation of the a-rotation
`v^{(a)}_j = pi((a+j) mod n)`). Let `S_a(b) = sigma_a(0) + ... + sigma_a(b-1)`
(`S_a(0) = 0`).

**C37.** For every `n >= 1`, every permutation `pi` of `Z_n`, and every
`a in Z_n`:

```
min_{b in Z_n} inv(a,b) = inv_count(sigma_a) - max_{0<=b<=n} [2*S_a(b) - b*(n-1)].
```

Consequently `I(pi) = min_a` of the right-hand side, computable in `O(n^2)`
time per `pi` (against `O(n^4)` for the brute double loop over `(a,b)` with a
naive `O(n^2)` inversion count), using the incremental identity

```
inv_count(sigma_{a+1}) = inv_count(sigma_a) + (n-1) - 2*pi(a)    (C37, step 2)
```

to get all `n` values `inv_count(sigma_a)` in `O(n)` after one `O(n log n)`
base computation, and an `O(n)` prefix-sum sweep of `sigma_a` per `a` for the
`max` term.

## 1. Proof

### Step 1 — telescoping in `b`

Fix `a` and write `v = v^{(a)}`, `w^{(b)} = w^{a,b}` (`w^{(b)}_j = (v_j-b) mod
n`). Compare `w^{(b)}` and `w^{(b+1)}`: every entry decreases by 1 mod n. Let
`j0` be the unique index with `v_{j0} = b`, i.e. `w^{(b)}_{j0} = 0`. For
`j != j0`, `w^{(b)}_j >= 1`, so `w^{(b+1)}_j = w^{(b)}_j - 1` (no wraparound);
for `j0`, `w^{(b+1)}_{j0} = n-1` (wraps from 0). So the sequence
`w^{(b+1)}` is obtained from `w^{(b)}` by decreasing every entry except
position `j0` by 1, and setting position `j0` to the new maximum `n-1`.

Relative order among indices `!= j0` is unchanged by a uniform decrement, so
none of those `C(n-1,2)` pairs change inversion status. Only pairs involving
`j0` can change:

- In `w^{(b)}`, since `w^{(b)}_{j0}=0` is the minimum, `j0` contributes an
  inversion with exactly the indices `j < j0` (i.e. `j0` "before" a larger
  value never inverts by being smallest at the front is not automatic —
  concretely: pair `(j,j0)` with `j<j0` is inverted iff `w_j > w_{j0} = 0`,
  true for all such `j` since all other entries are `>=1`; pair `(j0,k)` with
  `k>j0` is never inverted since `w_{j0}=0` is smallest). So `j0`
  contributes exactly `j0` inversions (one per index before it).
- In `w^{(b+1)}`, `w^{(b+1)}_{j0} = n-1` is the maximum, so symmetric
  reasoning gives `j0` contributes exactly `n-1-j0` inversions (one per
  index after it).

Hence

```
inv(a,b+1) - inv(a,b) = (n-1-j0) - j0 = (n-1) - 2*j0,    j0 = pos of value b in v^{(a)}.
```

By definition `v_{j0} = b` means `pi(a+j0 mod n) = b`, i.e.
`j0 = pinv(b) - a mod n = sigma_a(b)`. So

```
inv(a, b+1) - inv(a, b) = (n-1) - 2*sigma_a(b),   for b = 0..n-1.       (*)
```

### Step 2 — telescoping sum and the max

Summing (*) from `0` to `b-1`:

```
inv(a,b) = inv(a,0) + b*(n-1) - 2*S_a(b),   S_a(b) = sum_{t<b} sigma_a(t).
```

`inv(a,0) = inv_count(v^{(a)})`. A permutation and its inverse have equal
inversion counts (standard fact: `(i,j)` is an inversion of `v` iff
`(v_j,v_i)` is an inversion of `v^{-1}`, a bijection between the two
inversion sets), and `sigma_a` is exactly `(v^{(a)})^{-1}` (by the same
computation as above: `sigma_a(t) = pinv(t)-a mod n` is the index `j` with
`v^{(a)}_j = t`). So `inv(a,0) = inv_count(sigma_a)`.

Substituting and rearranging: `inv(a,b) = inv_count(sigma_a) - [2 S_a(b) -
b(n-1)]`. Minimizing over `b in {0,...,n-1}` (equivalently `b in {0,...,n}`,
since `S_a(n) = sum sigma_a = n(n-1)/2` gives `2 S_a(n) - n(n-1) = 0 =
2S_a(0) - 0`, so `b=n` and `b=0` tie) is the same as maximizing
`2 S_a(b) - b(n-1)` over that range, giving C37. ∎

### Step 3 — the `a`-increment identity

`sigma_{a+1}(t) = sigma_a(t) - 1 mod n` for every `t`, i.e. `sigma_{a+1}` is
the value-decrement of `sigma_a` by 1 — exactly the operation analyzed in
Step 1, with `sigma_a` playing the role of `v` and the roles of `n-1-j0`/`j0`
symmetric. The index `t0` with `sigma_a(t0) = 0` is (by the definition of
`sigma_a`) `t0 = pi(a)` (since `sigma_a(t)=0` means `pinv(t)=a`, i.e.
`t = pi(a)`). By the same computation as (*):

```
inv_count(sigma_{a+1}) - inv_count(sigma_a) = (n-1) - 2*pi(a).
```

∎ (matches C37 "step 2" above; this is what makes the whole `a`-sweep
`O(n)` after the first `O(n log n)` inversion count).

## 2. Scope and what this does *not* prove

C37 is an **identity**, true for every `pi` and every `a` (not a bound).
It reduces the search over `(a,b) in Z_n^2` for `I(pi)` to a search over
`a in Z_n` only, with an `O(n)` closed form for the inner `b`-minimum — an
algorithmic speedup (`O(n^2)` instead of `O(n^4)` per `pi`), not a proof of
H13-I. In particular:

- The stronger, `a`-free claim "`inv_count(sigma) - max_b[2S(b)-b(n-1)] <=
  floor((n-1)^2/4)` for every permutation `sigma`" is **false**: the Lemma C
  counterexample of `docs/notes/h13_line_model.md` §4 (`n=7`,
  `w=(0,5,4,3,2,1,6)`, `inv(w)=10`) has `min_b inv_count(rotate_b(w)) = 10 >
  9 = floor(36/4)` (checked directly), so the `b`-only optimization is not
  enough on its own — the search over `a` is still essential, consistent
  with the existing negative results on H13-I (`docs/notes/h13_line_model.md`
  §1.7, §4, and this session's `check_h13i_search.py`).
- C37 does not by itself bound `I(pi)`; it is a tool. `experiments/h13i_fast.py`
  uses it to extend the numerical ladder for H13-I past the `O(n^4)`-search
  range (was 4 <= n <= 10 for C33; this session's run reaches n = 11
  exhaustively and larger n by sampling, see `data/runs/h13i_fast/`).

## 3. Verification

`checks/check_C37.py`: (a) exact-vs-brute cross-check of the identity on
random permutations, `4 <= n <= 9` (500 trials, brute double loop over all
`(a,b)` with a naive `O(n^2)` inversion count); (b) `fast_I` (C37-based)
reproduces the brute `I(pi)` exactly on all permutations, `4 <= n <= 9`
(exhaustive); (c) `fast_I` reproduces the known `max_pi I(pi) =
floor((n-1)^2/4)` values of C33/C35 for `4 <= n <= 10` exhaustively, and
extends to `n = 11` (`data/runs/h13i_fast/report.json`). PASS.
