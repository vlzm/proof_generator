# circle_sort_equivalence report

I(pi) == noas(pi) for every pi, exhaustively, for n in [3, 7] (0 mismatches at every n); the easy direction I(pi) >= noas(pi) is additionally a proved inequality (assert never triggered). General n: open.

| n | checked | max_I | max_noas | floor((n-1)^2/4) | mismatches |
|---|---|---|---|---|---|
| 3 | 6 | 1 | 1 | 1 | 0 |
| 4 | 24 | 2 | 2 | 2 | 0 |
| 5 | 120 | 4 | 4 | 4 | 0 |
| 6 | 720 | 6 | 6 | 6 | 0 |
| 7 | 5040 | 9 | 9 | 9 | 0 |
