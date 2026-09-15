# check_C37/C38 -- arc/quadrant reformulation, single-shift rules ruled out

Команда: `python3 checks/check_C37.py --pairs_amax 9 --rules_amax 8`. Версии: check_C37-1.0, oracle-1.0.

```text
== check_C37-1.0 oracle-1.0 args={'pairs_amax': 9, 'rules_amax': 8}
C37 n=2: 2 pi x 4 (a,b) pairs x C(2,2) inner pairs -- exhaustive, 0.00 s
C37 n=3: 6 pi x 9 (a,b) pairs x C(3,2) inner pairs -- exhaustive, 0.00 s
C37 n=4: 24 pi x 16 (a,b) pairs x C(4,2) inner pairs -- exhaustive, 0.00 s
C37 n=5: 120 pi x 25 (a,b) pairs x C(5,2) inner pairs -- exhaustive, 0.02 s
C37 n=6: 720 pi x 36 (a,b) pairs x C(6,2) inner pairs -- exhaustive, 0.17 s
C37 n=7: 20 pi x 49 (a,b) pairs x C(7,2) inner pairs -- sampled, 0.01 s
C37 n=8: 20 pi x 64 (a,b) pairs x C(8,2) inner pairs -- sampled, 0.01 s
C37 n=9: 20 pi x 81 (a,b) pairs x C(9,2) inner pairs -- sampled, 0.02 s
C38 n=4 rule a=0: worst I_value_only = 2 vs bound 2 at (0, 2, 3, 1) -> within bound
C38 n=4 rule a=pos(0): worst I_value_only = 2 vs bound 2 at (0, 2, 3, 1) -> within bound
C38 n=4 rule a=pos(0)+1: worst I_value_only = 2 vs bound 2 at (0, 3, 1, 2) -> within bound
C38 n=4 rule a=pos(median): worst I_value_only = 2 vs bound 2 at (0, 2, 1, 3) -> within bound
C38 n=5 rule a=0: worst I_value_only = 4 vs bound 4 at (0, 2, 4, 3, 1) -> within bound
C38 n=5 rule a=pos(0): worst I_value_only = 4 vs bound 4 at (0, 2, 4, 3, 1) -> within bound
C38 n=5 rule a=pos(0)+1: worst I_value_only = 4 vs bound 4 at (0, 4, 1, 3, 2) -> within bound
C38 n=5 rule a=pos(median): worst I_value_only = 4 vs bound 4 at (0, 3, 2, 4, 1) -> within bound
C38 n=6 rule a=0: worst I_value_only = 6 vs bound 6 at (0, 2, 4, 5, 3, 1) -> within bound
C38 n=6 rule a=pos(0): worst I_value_only = 6 vs bound 6 at (0, 2, 4, 5, 3, 1) -> within bound
C38 n=6 rule a=pos(0)+1: worst I_value_only = 6 vs bound 6 at (0, 1, 5, 4, 3, 2) -> within bound
C38 n=6 rule a=pos(median): worst I_value_only = 6 vs bound 6 at (0, 1, 5, 3, 2, 4) -> within bound
C38 n=7 rule a=0: worst I_value_only = 10 vs bound 9 at (0, 5, 4, 3, 2, 1, 6) -> EXCEEDS (rule insufficient)
C38 n=7 rule a=pos(0): worst I_value_only = 10 vs bound 9 at (0, 5, 4, 3, 2, 1, 6) -> EXCEEDS (rule insufficient)
C38 n=7 rule a=pos(0)+1: worst I_value_only = 10 vs bound 9 at (0, 1, 6, 5, 4, 3, 2) -> EXCEEDS (rule insufficient)
C38 n=7 rule a=pos(median): worst I_value_only = 10 vs bound 9 at (0, 6, 5, 4, 2, 3, 1) -> EXCEEDS (rule insufficient)
C38 n=8 rule a=0: worst I_value_only = 13 vs bound 12 at (0, 6, 5, 3, 4, 2, 1, 7) -> EXCEEDS (rule insufficient)
C38 n=8 rule a=pos(0): worst I_value_only = 13 vs bound 12 at (0, 6, 5, 3, 4, 2, 1, 7) -> EXCEEDS (rule insufficient)
C38 n=8 rule a=pos(0)+1: worst I_value_only = 13 vs bound 12 at (0, 1, 7, 6, 4, 5, 3, 2) -> EXCEEDS (rule insufficient)
C38 n=8 rule a=pos(median): worst I_value_only = 13 vs bound 12 at (0, 6, 7, 5, 4, 2, 3, 1) -> EXCEEDS (rule insufficient)
C38 summary: rules with a confirmed counterexample (n <= 8): ['a=0', 'a=pos(0)', 'a=pos(0)+1', 'a=pos(median)']
verdict (C37 formula correctness): PASS
```
