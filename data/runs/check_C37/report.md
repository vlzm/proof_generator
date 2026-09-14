# check_C37 — замкнутая форма f(a, b) для двойного разреза (H13-I)

Команда: `python3 checks/check_C37.py --trials 200 --nmax 9`. Версия: check_C37-1.0.

```text
== check_C37-1.0 args={'trials': 200, 'nmax': 9, 'seed': 1}
part A+B: 7128 (pi, a, b) triples over 200 permutations (n <= 9) -> all match
part C n=4: max_pi I(pi) = 2 (bound floor((n-1)^2/4) = 2) -> ok, 0.00 s
  cross-check vs data/runs/line_profile/profile_n4.json: max_I=2 -> match
part C n=5: max_pi I(pi) = 4 (bound floor((n-1)^2/4) = 4) -> ok, 0.00 s
  cross-check vs data/runs/line_profile/profile_n5.json: max_I=4 -> match
part C n=6: max_pi I(pi) = 6 (bound floor((n-1)^2/4) = 6) -> ok, 0.01 s
  cross-check vs data/runs/line_profile/profile_n6.json: max_I=6 -> match
part C n=7: max_pi I(pi) = 9 (bound floor((n-1)^2/4) = 9) -> ok, 0.06 s
  cross-check vs data/runs/line_profile/profile_n7.json: max_I=9 -> match
part C n=8: max_pi I(pi) = 12 (bound floor((n-1)^2/4) = 12) -> ok, 0.59 s
  cross-check vs data/runs/line_profile/profile_n8.json: max_I=12 -> match
part C n=9: max_pi I(pi) = 16 (bound floor((n-1)^2/4) = 16) -> ok, 6.23 s
  cross-check vs data/runs/line_profile/profile_n9.json: max_I=16 -> match
verdict: PASS
```
