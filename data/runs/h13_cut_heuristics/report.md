# H13-I cut-heuristics (session 9)

Two single-feature candidates for choosing the position rotation a (then TRUE min over value shift b), tested exhaustively against target = floor((n-1)^2/4). worst_gap > 0 means the candidate fails to reach the target on the given (minimal, for that candidate) example.

| n | target | argmin worst_gap | argmin example | argmax worst_gap | argmax example | checked | s |
|---|---|---|---|---|---|---|---|
| 4 | 2 | 0 | (0, 3, 2, 1) | 0 | (0, 3, 1, 2) | 24 | 0.0 |
| 5 | 4 | 0 | (0, 4, 3, 2, 1) | 0 | (0, 4, 1, 3, 2) | 120 | 0.0 |
| 6 | 6 | 0 | (0, 4, 3, 2, 1, 5) | 0 | (0, 5, 1, 3, 4, 2) | 720 | 0.02 |
| 7 | 9 | 1 | (0, 5, 4, 3, 2, 1, 6) | 0 | (0, 6, 1, 3, 5, 4, 2) | 5040 | 0.12 |
| 8 | 12 | 0 | (0, 5, 4, 3, 2, 7, 1, 6) | 0 | (0, 7, 1, 3, 5, 6, 4, 2) | 40320 | 0.92 |
| 9 | 16 | 1 | (0, 6, 5, 4, 3, 2, 8, 1, 7) | 1 | (0, 8, 1, 6, 5, 4, 3, 2, 7) | 362880 | 8.03 |
| 10 | 20 | 2 | (0, 7, 6, 5, 4, 3, 2, 1, 9, 8) | 1 | (0, 9, 1, 7, 6, 4, 5, 3, 2, 8) | 3628800 | 102.87 |
