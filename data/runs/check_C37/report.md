# check_C37 — PASS

Версия check_C37-1.0, seed 1, 8.4 с.

```
check_C37-1.0  seed=1
part A: lemmas 0-5 (C37)
  A n=4 exhaustive 24 permutations, 0.00 s
  A n=5 exhaustive 120 permutations, 0.03 s
  A n=6 exhaustive 720 permutations, 0.31 s
  A n=7 exhaustive 5040 permutations, 3.27 s
  A n=8 sampled 200 permutations, 0.20 s
  A n=9 sampled 200 permutations, 0.28 s
  A n=10 sampled 200 permutations, 0.39 s
  A n=11 sampled 200 permutations, 0.52 s
part B: lemma 6, reflection grid (C37)
  B reflections checked for 4 <= n <= 24
part C: C38, sum_b min_a F(a,b) <= n*M
  C n=4 reps=6 max sum=8 (bound 8) argmax=[[0, 3, 2, 1]] maxI=2
  C n=5 reps=24 max sum=20 (bound 20) argmax=[[0, 4, 3, 2, 1]] maxI=4
  C n=6 reps=120 max sum=36 (bound 36) argmax=[[0, 5, 4, 3, 2, 1]] maxI=6
  C n=7 reps=720 max sum=63 (bound 63) argmax=[[0, 6, 5, 4, 3, 2, 1]] maxI=9
  C n=8 reps=5040 max sum=96 (bound 96) argmax=[[0, 7, 6, 5, 4, 3, 2, 1]] maxI=12
  C external scan rows read: n = [4, 5, 6, 7, 8, 9, 10, 11, 12]
part D: C39, circular swap distance vs I
  D n=4 max cs=2 (M=2) #(cs<I)=0 first=None
  D n=5 max cs=4 (M=4) #(cs<I)=0 first=None
  D n=6 max cs=6 (M=6) #(cs<I)=0 first=None
  D n=7 max cs=9 (M=9) #(cs<I)=0 first=None
  D n=8 max cs=12 (M=12) #(cs<I)=16 first=([0, 3, 6, 1, 4, 7, 2, 5], 8, 9)
part E: smallest counterexamples of the refuted candidates
  E smallest counterexamples reproduced: ['T1', 'T2', 'T3', 'T4', 'anti', 'diag', 'linear', 'rowmax']
```
