# H13-I candidate strategies -- report

Exhaustive over all pi, n = 4..9. Version h13i_candidates-1.0.

| n | bound floor((n-1)^2/4) | max I(pi) | max minq_avgs (vs n*bound) | max anchor_min (vs bound) | min best_align x (need ceil(n/2)) | perfect alignments / perms | time |
|---|---|---|---|---|---|---|---|
| 4 | 2 | 2 | 14 (8) | 6 (2) | 2 (2) | 24/24 | 0.0s |
| 5 | 4 | 4 | 30 (20) | 10 (4) | 2 (3) | 110/120 | 0.0s |
| 6 | 6 | 6 | 55 (36) | 15 (6) | 2 (3) | 516/720 | 0.1s |
| 7 | 9 | 9 | 91 (63) | 21 (9) | 3 (4) | 3346/5040 | 1.2s |
| 8 | 12 | 12 | 140 (96) | 28 (12) | 3 (4) | 15856/40320 | 13.5s |
| 9 | 16 | 16 | 204 (144) | 36 (16) | 3 (5) | 133560/362880 | 177.7s |

Worst cases (per n):

- n=4: minq_avgs worst pi = [0, 3, 2, 1]; anchor worst pi = [0, 3, 2, 1]; align worst pi = [0, 1, 2, 3]
- n=5: minq_avgs worst pi = [0, 4, 3, 2, 1]; anchor worst pi = [0, 4, 3, 2, 1]; align worst pi = [0, 2, 4, 1, 3]
- n=6: minq_avgs worst pi = [0, 5, 4, 3, 2, 1]; anchor worst pi = [0, 5, 4, 3, 2, 1]; align worst pi = [0, 1, 3, 5, 2, 4]
- n=7: minq_avgs worst pi = [0, 6, 5, 4, 3, 2, 1]; anchor worst pi = [0, 6, 5, 4, 3, 2, 1]; align worst pi = [0, 1, 3, 4, 6, 2, 5]
- n=8: minq_avgs worst pi = [0, 7, 6, 5, 4, 3, 2, 1]; anchor worst pi = [0, 7, 6, 5, 4, 3, 2, 1]; align worst pi = [0, 1, 2, 4, 5, 7, 3, 6]
- n=9: minq_avgs worst pi = [0, 8, 7, 6, 5, 4, 3, 2, 1]; anchor worst pi = [0, 8, 7, 6, 5, 4, 3, 2, 1]; align worst pi = [0, 2, 4, 6, 8, 1, 3, 5, 7]
