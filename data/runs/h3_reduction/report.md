# h3_reduction — индукция H3: бюджет n−1 против одноэлементной редукции

Команда: `python3 experiments/h3_reduction.py --nmax 8 --sample-from 9 --sample 100000 --prep-max 8 --seed 20260915`. Версии: h3_reduction-1.0, oracle-1.0, таблицы C2v.

```text
== h3_reduction-1.0 core=oracle-1.0 args={'nmax': 8, 'sample_from': 9, 'sample': 100000, 'prep_max': 8, 'seed': 20260915}
n=5 (exhaustive, 120 states, tables C2v): max [d_n(p) - min_k d_(n-1)(reduce(p,k))] = 5 vs n-1 = 4 -> > n-1: the cheapest reduction is incompatible with a 1:1 simulation; worst p = (0, 1, 2, 4, 3); 0.0 s
   histogram of the excess: 0:1, 1:9, 2:17, 3:20, 4:44, 5:29
   best fixed deletion rule (delete array index k): k=0:5, k=1:5, k=2:5, k=3:5, k=4:5  (budget 4)
   slack (zero preparation): min over states of [max_k d_(n-1) - (d_n - (n-1))] = 0 at p = (0, 1, 4, 3, 2); states with negative slack: 0 of 120
   preparation (all words |u| <= 4): rescued 0 of 0 negative-slack states; unrescued: 0
n=6 (exhaustive, 720 states, tables C2v): max [d_n(p) - min_k d_(n-1)(reduce(p,k))] = 8 vs n-1 = 5 -> > n-1: the cheapest reduction is incompatible with a 1:1 simulation; worst p = (0, 1, 2, 4, 5, 3); 0.0 s
   histogram of the excess: 0:1, 1:10, 2:28, 3:62, 4:131, 5:238, 6:204, 7:30, 8:16
   best fixed deletion rule (delete array index k): k=0:7, k=1:7, k=2:8, k=3:8, k=4:8, k=5:8  (budget 5)
   slack (zero preparation): min over states of [max_k d_(n-1) - (d_n - (n-1))] = 0 at p = (0, 1, 4, 5, 2, 3); states with negative slack: 0 of 720
   preparation (all words |u| <= 5): rescued 0 of 0 negative-slack states; unrescued: 0
n=7 (exhaustive, 5040 states, tables C2v): max [d_n(p) - min_k d_(n-1)(reduce(p,k))] = 9 vs n-1 = 6 -> > n-1: the cheapest reduction is incompatible with a 1:1 simulation; worst p = (0, 1, 2, 4, 5, 3, 6); 0.1 s
   histogram of the excess: 0:1, 1:13, 2:72, 3:84, 4:514, 5:1088, 6:1639, 7:1166, 8:459, 9:4
   best fixed deletion rule (delete array index k): k=0:8, k=1:8, k=2:8, k=3:9, k=4:8, k=5:9, k=6:8  (budget 6)
   slack (zero preparation): min over states of [max_k d_(n-1) - (d_n - (n-1))] = -1 at p = (6, 1, 4, 5, 2, 3, 0); states with negative slack: 1 of 5040
   preparation (all words |u| <= 6): rescued 1 of 1 negative-slack states; unrescued: 0
n=8 (exhaustive, 40320 states, tables C2v): max [d_n(p) - min_k d_(n-1)(reduce(p,k))] = 11 vs n-1 = 7 -> > n-1: the cheapest reduction is incompatible with a 1:1 simulation; worst p = (0, 1, 2, 3, 5, 6, 7, 4); 0.9 s
   histogram of the excess: 0:1, 1:14, 2:94, 3:247, 4:919, 5:3052, 6:8987, 7:10946, 8:9792, 9:5732, 10:492, 11:44
   best fixed deletion rule (delete array index k): k=0:10, k=1:10, k=2:11, k=3:11, k=4:11, k=5:11, k=6:11, k=7:11  (budget 7)
   slack (zero preparation): min over states of [max_k d_(n-1) - (d_n - (n-1))] = -2 at p = (2, 3, 7, 4, 0, 1, 5, 6); states with negative slack: 62 of 40320
   preparation (all words |u| <= 7): rescued 62 of 62 negative-slack states; unrescued: 0
n=9 (random sample of 100000, 100000 states, tables C2v): max [d_n(p) - min_k d_(n-1)(reduce(p,k))] = 12 vs n-1 = 8 -> > n-1: the cheapest reduction is incompatible with a 1:1 simulation; worst p = (1, 2, 3, 5, 6, 7, 4, 8, 0); 3.4 s
   histogram of the excess: 0:1, 1:5, 2:67, 3:120, 4:630, 5:3577, 6:7695, 7:20081, 8:25791, 9:26961, 10:11562, 11:3507, 12:3
   best fixed deletion rule (delete array index k): k=0:11, k=1:11, k=2:11, k=3:11, k=4:11, k=5:11, k=6:12, k=7:12, k=8:11  (budget 8)
   slack (zero preparation): min over states of [max_k d_(n-1) - (d_n - (n-1))] = -2 at p = (8, 1, 4, 6, 7, 2, 3, 5, 0); states with negative slack: 19 of 100000
verdict: OBSTRUCTION: states with negative slack exist
```
