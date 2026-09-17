# check_C37 — рекуррента для инверсий одного вращения значений (H13-I, сессия 9)

Команда: `python3 checks/check_C37.py --emax 8 --rand_nmax 200 --rand_trials 300 --seed 0`. Версии: check_C37-1.0, oracle-1.0.

```text
== check_C37-1.0 oracle-1.0 args={'emax': 8, 'rand_nmax': 200, 'rand_trials': 300, 'seed': 0}
recursion exhaustive n=2: 2 permutations, 0.00 s
recursion exhaustive n=3: 6 permutations, 0.00 s
recursion exhaustive n=4: 24 permutations, 0.00 s
recursion exhaustive n=5: 120 permutations, 0.00 s
recursion exhaustive n=6: 720 permutations, 0.02 s
recursion exhaustive n=7: 5040 permutations, 0.21 s
recursion exhaustive n=8: 40320 permutations, 2.19 s
recursion random n=9: 300 trials, 0.02 s
recursion random n=10: 300 trials, 0.03 s
recursion random n=20: 300 trials, 0.11 s
recursion random n=50: 300 trials, 0.73 s
recursion random n=200: 300 trials, 13.61 s
inv(v)=inv(v^-1) n=2: 2 permutations checked
inv(v)=inv(v^-1) n=3: 6 permutations checked
inv(v)=inv(v^-1) n=4: 24 permutations checked
inv(v)=inv(v^-1) n=5: 120 permutations checked
inv(v)=inv(v^-1) n=6: 720 permutations checked
inv(v)=inv(v^-1) n=7: 5040 permutations checked
inv(v)=inv(v^-1) n=8: 40320 permutations checked
value-rotation-only n=7 v=(0, 5, 4, 3, 2, 1, 6): min_b inv = 10, floor((n-1)^2/4) = 9, exceeds by 1 -> ok (matches recorded failure)
value-rotation-only n=8 v=(0, 6, 5, 3, 4, 2, 1, 7): min_b inv = 13, floor((n-1)^2/4) = 12, exceeds by 1 -> ok (matches recorded failure)
verdict: PASS
```
