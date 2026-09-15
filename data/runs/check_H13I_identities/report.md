# check_H13I_identities report (check_H13I_identities-1.0, oracle-1.0)

Overall: PASS (total_fails=0)

## Part A/B: shift identities (C37), duality (C37), double counting (C38)

| n | mode | count | c_shift fails | q_shift fails | duality fails | double_counting fails | seconds |
|---|---|---|---|---|---|---|---|
| 4 | exhaustive | 24 | 0 | 0 | 0 | 0 | 0.0 |
| 5 | exhaustive | 120 | 0 | 0 | 0 | 0 | 0.0 |
| 6 | exhaustive | 720 | 0 | 0 | 0 | 0 | 0.4 |
| 7 | exhaustive | 5040 | 0 | 0 | 0 | 0 | 3.8 |
| 8 | sampled (60 random + sigma_n + rev_n) | 62 | 0 | 0 | 0 | 0 | 0.1 |
| 9 | sampled (60 random + sigma_n + rev_n) | 62 | 0 | 0 | 0 | 0 | 0.1 |
| 10 | sampled (60 random + sigma_n + rev_n) | 62 | 0 | 0 | 0 | 0 | 0.2 |
| 11 | sampled (60 random + sigma_n + rev_n) | 62 | 0 | 0 | 0 | 0 | 0.2 |
| 12 | sampled (60 random + sigma_n + rev_n) | 62 | 0 | 0 | 0 | 0 | 0.3 |
| 13 | sampled (60 random + sigma_n + rev_n) | 62 | 0 | 0 | 0 | 0 | 0.4 |
| 14 | sampled (60 random + sigma_n + rev_n) | 62 | 0 | 0 | 0 | 0 | 0.5 |
| 15 | sampled (60 random + sigma_n + rev_n) | 62 | 0 | 0 | 0 | 0 | 0.6 |
| 16 | sampled (60 random + sigma_n + rev_n) | 62 | 0 | 0 | 0 | 0 | 0.7 |
| 17 | sampled (60 random + sigma_n + rev_n) | 62 | 0 | 0 | 0 | 0 | 0.9 |
| 18 | sampled (60 random + sigma_n + rev_n) | 62 | 0 | 0 | 0 | 0 | 1.0 |

## Part C: why averaging and the c=0 family fail (record, no new claim)

| n | bound | avg inv sigma_n | avg inv rev_n | averaging fails? | max_pi min_q inv(.,q,0) | c=0 family sufficient? |
|---|---|---|---|---|---|---|
| 4 | 2 | 7/2 | 7/2 | True | 4 | False |
| 5 | 4 | 6 | 6 | True | 6 | False |
| 6 | 6 | 55/6 | 55/6 | True | 10 | False |
| 7 | 9 | 13 | 13 | True | 12 | False |
| 8 | 12 | 35/2 | 35/2 | True | 17 | False |
| 9 | 16 | 68/3 | 68/3 | True | None | None |
| 10 | 20 | 57/2 | 57/2 | True | None | None |
