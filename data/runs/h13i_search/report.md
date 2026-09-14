# H13-I search — report (session 9, 14.09.2026)

Goal: prove H13-I (`I(pi) <= floor((n-1)^2/4)` for all n >= 4, pi in S_n;
PLAN.md H13, `docs/notes/h13_line_model.md` S6), budget 1 session until
written verdict. Verdict: not proved, not refuted; three candidate
approaches rejected with counterexamples; numerical range extended.
Full writeup: `docs/notes/h13i_verdict.md`.

## Command / versions

- `gcc -O2 -o h13i_toric experiments/h13i_toric.c` (h13i_toric-1.0),
  core oracle-1.0 (no distance table needed — H13-I is purely about I(pi)).
- `python3 experiments/h13i_candidates.py --nmax 8` (h13i_candidates-1.0).
- `python3 experiments/h13i_sample.py --nlist 11,12,15,20,25,30,40,50,60,80,100,150,200 --seed 42 --random-trials 8` (h13i_sample-1.0).

## Coverage and results

1. **Exhaustive**, all permutations, all n^2 cuts: 4 <= n <= 11 (extends C33's
   prior 4 <= n <= 10). `toric_n4_11.log`: `max_I(pi) == floor((n-1)^2/4)`
   exactly at every n, achieved by a reflection (the lexicographically-first
   is `(0, n-1, n-2, ..., 1)`); no counterexample. n = 11: 39 916 800
   permutations, 2 m 15 s on one core.
2. **Exhaustive n = 12**: 4-way sharded (`skip_mod=4`), ~9 min wall-clock on
   4 cores (matches the ~14x scaling estimate from n=10 -> n=11). All
   4 x 119 750 400 = 479 001 600 = 12! permutations checked (sum of the
   `checked` fields across shards equals 12! exactly); `max_I = 30 =
   floor(11^2/4)`, no counterexample. Logs: `toric_n12_shard{0,1,2,3}.log`.
3. **SAMPLED n = 11..200** (structured families + 8 random permutations per
   n): `sample_report.md`/`.json`. No counterexample; reflections hit the
   target exactly at every n, other families stay strictly below with a
   growing margin.
4. **Rejected candidate cut-selection rules** (would have given an O(n),
   not O(n^2), certificate): `candidates_report.md`/`.json`. Both the
   "diagonal" rule (cut origin = a point of pi itself) and the
   "mode-shift" rule (fix c to the corrected mode of `pi(i)-i`, vary only q)
   fail starting at n = 4, exhaustively confirmed 3 <= n <= 8.
5. **Full-grid averaging re-derived and confirmed to fail** (already noted
   in session 8; here with an explicit exact number): average of inv(q,c)
   over all n^2 cuts on the n = 10 reflection is 28.5, strictly above the
   target 20 — so `min <= average` cannot certify H13-I in general.

## Conclusion

No new counterexample at any tested n (up to 200, exhaustively to n = 12);
reflections remain the unique known extremal family, consistent with C33.
Three concrete proof strategies are now recorded as refuted (not merely
"not attempted"). Next step per the verdict note: induction H3 (PLAN S4.3
p.3), not further search for an O(n)/O(n^2) certificate rule along these
particular lines.
