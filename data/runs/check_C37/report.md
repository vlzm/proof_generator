# check_C37 — отчёт

check_C37-1.0, hmax=8, 22.3 с.

```text
L1: equivalence inv <= floor((n-1)^2/4) <=> S >= floor(n/2), 2<=n<=200: PASS
L2: cut family = switchings by A_q xor B_s, one step = one vertex, all pi and cuts 4<=n<=7: PASS
T: n=4, all 64 signings: min over signings of max over switchings S = 2 (need >= 2)
T: n=5, all 1024 signings: min over signings of max over switchings S = 2 (need >= 2)
T: n=6, permutation signings: min max_switch S = 3 (need >= 3); hist of (max_switch - max_cut) = {0: 120}
T: n=7, permutation signings: min max_switch S = 3 (need >= 3); hist of (max_switch - max_cut) = {0: 720}
T: n=8, permutation signings: min max_switch S = 4 (need >= 4); hist of (max_switch - max_cut) = {0: 5014, 2: 25, 6: 1}
T: all-minus class maximum = floor(n/2) exactly, 2<=n<=40; theorem: PASS
L7: identity S = (4/n)K + (2/n^2)T - 1 - n/2, all pi and cuts 4<=n<=8: PASS
L7: mean number of non-inversions over all n^2 cuts of a reflection = (n^2-1)/6 < floor(n^2/4), 4<=n<=12: PASS
P6: n=4, representatives with no stable cut: 0 
P6: n=5, representatives with no stable cut: 0 
P6: n=6, representatives with no stable cut: 0 
P6: n=7, representatives with no stable cut: 0 
P6: n=8, representatives with no stable cut: 2 [[0, 3, 6, 1, 4, 7, 2, 5], [0, 5, 2, 7, 4, 1, 6, 3]]
P6: no stable cut for 4<=n<=7: none; n=8: exactly pi(i)=3i and 5i: PASS
V: n=5: min over pi of max over (a,b) of sum_j S(j+a, pi(j)+b) = 8 < n*floor(n/2) = 10 at pi=[0, 1, 4, 3, 2]
V: n=7: min over pi of max over (a,b) of sum_j S(j+a, pi(j)+b) = 13 < n*floor(n/2) = 21 at pi=[0, 2, 1, 6, 5, 4, 3]
R: reflections -- line form, inv = C(c+1,2)+C(n-1-c,2), min over cuts = floor((n-1)^2/4), all-minus signing present, 4<=n<=24, all h: PASS
H: n=4 reps=6 max_I=2 (bound 2) argmax=1 classes, exactly the reflections: True [0.0s]
H: n=5 reps=24 max_I=4 (bound 4) argmax=1 classes, exactly the reflections: True [0.0s]
H: n=6 reps=120 max_I=6 (bound 6) argmax=1 classes, exactly the reflections: True [0.0s]
H: n=7 reps=720 max_I=9 (bound 9) argmax=1 classes, exactly the reflections: True [0.1s]
H: n=8 reps=5040 max_I=12 (bound 12) argmax=1 classes, exactly the reflections: True [1.0s]
H: families (reflections, affine, random, block) n<=40: min S over family always >= floor(n/2): PASS
```
