# H13-I proof attempts: report

h13i_attempts-1.0. Exhaustive 4 <= n <= 9. All pi.

| n | target | fixed_q fails | fixed_q worst gap | example | w0 fails | w0 worst gap | example | local-search fails | local-search worst gap | example |
|---|---|---|---|---|---|---|---|---|---|---|
| 4 | 2 | 0 | 0 | None | 4 | 2 | (0, 2, 1, 3) | 2 | 1 | (1, 3, 0, 2) |
| 5 | 4 | 0 | 0 | None | 13 | 3 | (0, 3, 1, 4, 2) | 15 | 2 | (1, 3, 0, 4, 2) |
| 6 | 6 | 0 | 0 | None | 62 | 3 | (0, 3, 4, 1, 2, 5) | 101 | 4 | (2, 5, 4, 1, 0, 3) |
| 7 | 9 | 7 | 1 | (0, 5, 4, 3, 2, 1, 6) | 215 | 5 | (3, 0, 4, 1, 5, 2, 6) | 676 | 5 | (2, 5, 6, 4, 1, 0, 3) |
| 8 | 12 | 16 | 1 | (0, 6, 5, 3, 4, 2, 1, 7) | 920 | 4 | (3, 0, 4, 1, 5, 2, 6, 7) | 5240 | 7 | (2, 6, 7, 5, 1, 4, 0, 3) |
| 9 | 16 | 162 | 1 | (0, 2, 7, 6, 5, 4, 3, 8, 1) | 4109 | 6 | (4, 0, 5, 1, 6, 2, 7, 3, 8) | 47949 | 9 | (2, 6, 8, 7, 5, 1, 4, 0, 3) |

Conclusion: all three candidate reductions of H13-I fail exhaustively from small n (see docs/notes/h13_line_model.md §7 for the verdict). The full 2-parameter search over (q, c) remains necessary and no closed-form / single-pass construction was found this session.
