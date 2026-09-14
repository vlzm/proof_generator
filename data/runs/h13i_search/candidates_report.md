# H13-I candidate cut-selection rules -- report

Version h13i_candidates-1.0, core oracle-1.0. Exhaustive per n, 10.5 s total.

| n | target | worst_full | worst_diagonal | worst_mode_shift |
|---|---|---|---|---|
| 3 | 1 | 1 | 1 | 1 |
| 4 | 2 | 2 | 3 | 3 |
| 5 | 4 | 4 | 6 | 4 |
| 6 | 6 | 6 | 10 | 7 |
| 7 | 9 | 9 | 15 | 9 |
| 8 | 12 | 12 | 21 | 15 |

Both O(n) rules fail (worst > target) starting at n = 4; example permutations and their candidate vs. true I(pi) are in candidates_report.json under `examples`.
