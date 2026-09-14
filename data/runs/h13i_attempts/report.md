# H13-I proof attempts (h13i_attempts-1.0)

Four candidate simplifications, each tested exhaustively against `floor((n-1)^2/4)`. All four fail (counterexamples appear by n = 9 at the latest); see docs/notes/h13_line_model.md §7.

| n | bound | fix_a0 worst | cut_at_data worst | footrule worst | clean_interval worst |
|---|---|---|---|---|---|
| 4 | 2 | 2 | 3 | 4 | 2 |
| 5 | 4 | 4 | 6 | 6 | 6 |
| 6 | 6 | 6 | 10 | 8 | 10 |
| 7 | 9 | 10 | 15 | 12 | 15 |
| 8 | 12 | 13 | 21 | 16 | 21 |
| 9 | 16 | 17 | 28 | - | - |
