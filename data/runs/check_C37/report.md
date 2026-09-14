# check_C37-1.0 report

Overall: **PASS**  (elapsed 5.6 s)

## Part A -- cut-exchange recurrence

Checked 574488 (a,b) recurrence instances, coverage: [('exhaustive', 4), ('exhaustive', 5), ('exhaustive', 6), ('exhaustive', 7), ('random', 8, 30), ('random', 9, 30), ('random', 10, 30), ('random', 11, 30)]. Result: PASS

## Part B -- restricted one-parameter families (a=0 only; diagonal b=pi(a))

| n | target | a=0-only worst | a=0-only #exceeding | diagonal worst | diagonal #exceeding |
|---|---|---|---|---|---|
| 4 | 2 | 2 | 0/24 | 3 | 4/24 |
| 5 | 4 | 4 | 0/120 | 6 | 5/120 |
| 6 | 6 | 6 | 0/720 | 10 | 6/720 |
| 7 | 9 | 10 | 7/5040 | 15 | 56/5040 |

## Part C -- full 2D search sanity re-check (not a new claim; cross-check of C33)

| n | target | worst over all pi of min_{a,b} inv | OK |
|---|---|---|---|
| 4 | 2 | 2 | True |
| 5 | 4 | 4 | True |
| 6 | 6 | 6 | True |
| 7 | 9 | 9 | True |
