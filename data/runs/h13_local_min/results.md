# h13_local_min data (session 9)

Purpose: test whether 2D ("king move") local optimality of the double cut
(q, c), reparametrized as (q, S) with S the value-shift, at a fixed radius R
certifies `inv <= floor((n-1)^2/4)` (candidate path to H13-I). See
`docs/notes/h13_line_model.md` §7 for the write-up and `experiments/h13_local_min.c`
for the method (fast O(n^2)-per-permutation recurrence, checked against brute
force inversion counts for 4 <= n <= 7 via the `check` flag).

Command: `gcc -O2 -o experiments/h13_local_min experiments/h13_local_min.c`,
then `./experiments/h13_local_min n [check] R`. Raw stdout: `raw.txt`.

## Results (exhaustive over all n! permutations)

| n | bound = floor((n-1)^2/4) | R=1 worst local min | R=2 worst local min | R=3 worst local min |
|---|---|---|---|---|
| 4 | 2 | 2 | 2 | - |
| 5 | 4 | 4 | 4 | - |
| 6 | 6 | 7 | 6 | - |
| 7 | 9 | 11 | 9 | - |
| 8 | 12 | 16 | 12 | - |
| 9 | 16 | 22 | 16 | - |
| 10 | 20 | 29 | 21 | 20 |
| 11 | 25 | 37 | 29 | (not run) |

R = 1 (8 neighbours) already fails at n = 6; R = 2 (24 neighbours) matches the
bound exactly for 4 <= n <= 9 and fails at n = 10 and n = 11; R = 3 recovers
n = 10 but was not run at n = 11 (cost). No fixed radius tested here holds for
all n in the tested range; the radius needed, if bounded at all, is not
determined by this experiment.

Counterexample at n = 10, R = 2: `pi = (0, 1, 4, 9, 3, 8, 7, 6, 2, 5)`
(lexicographic rank 14638). Cell `(q, S) = (8, 4)` has `inv = 21`; no cell
within king-move distance 2 is smaller, but a cell at distance exactly 3
(e.g. `(4, 5)` or `(7, 0)`) has `inv = 14 = I(pi)`, the global minimum.

## Coverage

- Fast recurrence checked against brute-force inversion counting on all
  cells, all permutations, 4 <= n <= 7 (`check` flag; no mismatches).
- R = 1, 2 exhaustive over all n! permutations, 4 <= n <= 10 (n = 10: ~16-20s
  each); R = 1, 2 exhaustive at n = 11 (~2-4 min each depending on load).
- R = 3 exhaustive at n = 10 only (~20s); not run at n = 11 (would cost
  roughly 2x the R=2 n=11 time for a larger neighbourhood, decided not
  worth it without a specific hypothesis to test).
