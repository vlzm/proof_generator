# check_H13I -- алгоритм и конечные данные (сессия 9)

Команда: `python3 checks/check_H13I.py`. Версия check_H13I-1.0.

```text
== check_H13I-1.0 args={'n11': True, 'search_n': [12, 13, 14], 'search_iters': 1500, 'search_trials': 4}
n=4: line_toric_max_i (C, O(n^3)) max_I=2 floor((n-1)^2/4)=2 count=24 0.0s
n=4: brute force (independent Python, O(n^4)) max_I=2 vs fast C (O(n^3)) max_I=2 -> MATCH
n=5: line_toric_max_i (C, O(n^3)) max_I=4 floor((n-1)^2/4)=4 count=120 0.0s
n=5: brute force (independent Python, O(n^4)) max_I=4 vs fast C (O(n^3)) max_I=4 -> MATCH
n=6: line_toric_max_i (C, O(n^3)) max_I=6 floor((n-1)^2/4)=6 count=720 0.0s
n=6: brute force (independent Python, O(n^4)) max_I=6 vs fast C (O(n^3)) max_I=6 -> MATCH
n=7: line_toric_max_i (C, O(n^3)) max_I=9 floor((n-1)^2/4)=9 count=5040 0.0s
n=7: brute force (independent Python, O(n^4)) max_I=9 vs fast C (O(n^3)) max_I=9 -> MATCH
n=8: line_toric_max_i (C, O(n^3)) max_I=12 floor((n-1)^2/4)=12 count=40320 0.0s
n=8: brute force (independent Python, O(n^4)) max_I=12 vs fast C (O(n^3)) max_I=12 -> MATCH
n=9: line_toric_max_i (C, O(n^3)) max_I=16 floor((n-1)^2/4)=16 count=362880 0.4s
n=10: line_toric_max_i (C, O(n^3)) max_I=20 floor((n-1)^2/4)=20 count=3628800 5.5s
n=11: line_toric_max_i (C, O(n^3)) max_I=25 floor((n-1)^2/4)=25 count=39916800 73.3s
n=12: local search (annealing, targeted at maximizing I) best=30 target=30 exceeds_target=False 4.9s
n=13: local search (annealing, targeted at maximizing I) best=31 target=36 exceeds_target=False 6.4s
n=14: local search (annealing, targeted at maximizing I) best=42 target=42 exceeds_target=False 8.2s
verdict: PASS (this checker validates the algorithm and re-confirms no counterexample found in the ranges tested; H13-I itself remains CONJECTURED)
```
