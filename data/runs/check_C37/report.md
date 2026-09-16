# check_C37 — алгебра двойного разреза и граница свободного разреза

Команда: `python3 checks/check_C37.py --amax 8 --free-amax 8`. Версии: check_C37-1.0, oracle-1.0.

```text
== check_C37-1.0 oracle-1.0 args={'amax': 8, 'free_amax': 8}
lemma 1 n=4: all perms, (q,c) and (m,k) families coincide -> ok
lemma 1 n=5: all perms, (q,c) and (m,k) families coincide -> ok
lemma 1 n=6: all perms, (q,c) and (m,k) families coincide -> ok
lemma 2,3 n=4: all perms x 16 cuts -> ok (0.0 s)
lemma 2,3 n=5: all perms x 25 cuts -> ok (0.0 s)
lemma 2,3 n=6: all perms x 36 cuts -> ok (0.2 s)
lemma 2,3 n=7: all perms x 49 cuts -> ok (2.5 s)
lemma 2 n=8: 100 random (pi, m, k) -> ok
lemma 2 n=12: 100 random (pi, m, k) -> ok
lemma 2 n=20: 100 random (pi, m, k) -> ok
lemma 2 n=33: 100 random (pi, m, k) -> ok
lemma 2 n=50: 100 random (pi, m, k) -> ok
lemma 4 n=4: all perms, pair counts and Sigma(pi) -> ok
lemma 4 n=5: all perms, pair counts and Sigma(pi) -> ok
lemma 4 n=6: all perms, pair counts and Sigma(pi) -> ok
lemma 4: closed forms Sigma(id_n), Sigma(sigma_n) = Sigma(rev_n) and the averaging verdicts, 4 <= n <= 60 -> ok
theorem 5 n=4: all perms, min_pi max_x Q(x) = 2 >= floor(n/2) = 2 (cross-checked against the direct sum) -> ok (0.0 s)
theorem 5 n=5: all perms, min_pi max_x Q(x) = 2 >= floor(n/2) = 2 (cross-checked against the direct sum) -> ok (0.0 s)
theorem 5 n=6: all perms, min_pi max_x Q(x) = 3 >= floor(n/2) = 3 (cross-checked against the direct sum) -> ok (0.1 s)
theorem 5 n=7: all perms, min_pi max_x Q(x) = 3 >= floor(n/2) = 3 -> ok (0.7 s)
theorem 5 n=8: all perms, min_pi max_x Q(x) = 4 >= floor(n/2) = 4 -> ok (13.7 s)
theorem 5 n=4: 300 random +-1 signings of K_n, min max_x Q = 2 >= 2 -> ok
theorem 5 n=5: 300 random +-1 signings of K_n, min max_x Q = 2 >= 2 -> ok
theorem 5 n=6: 300 random +-1 signings of K_n, min max_x Q = 5 >= 3 -> ok
theorem 5 n=7: 300 random +-1 signings of K_n, min max_x Q = 7 >= 3 -> ok
theorem 5 n=8: 300 random +-1 signings of K_n, min max_x Q = 10 >= 4 -> ok
theorem 5 n=9: 300 random +-1 signings of K_n, min max_x Q = 10 >= 4 -> ok
theorem 5 n=10: 300 random +-1 signings of K_n, min max_x Q = 13 >= 5 -> ok
theorem 5: tightness eps = -1 (rev_n), max_x Q = floor(n/2), 4 <= n <= 14 -> ok
theorem 5 n=12: descent on 25 inputs, max resulting inv = 30 <= floor((n-1)^2/4) = 30 -> ok
theorem 5 n=20: descent on 27 inputs, max resulting inv = 90 <= floor((n-1)^2/4) = 90 -> ok
theorem 5 n=33: descent on 35 inputs, max resulting inv = 256 <= floor((n-1)^2/4) = 256 -> ok
theorem 5 n=50: descent on 35 inputs, max resulting inv = 600 <= floor((n-1)^2/4) = 600 -> ok
theorem 5 n=80: descent on 35 inputs, max resulting inv = 1560 <= floor((n-1)^2/4) = 1560 -> ok
theorem 5 n=101: descent on 35 inputs, max resulting inv = 2500 <= floor((n-1)^2/4) = 2500 -> ok
cor 5.2 n=4: implication holds; perms with no cut satisfying the criterion: 0
cor 5.2 n=5: implication holds; perms with no cut satisfying the criterion: 0
cor 5.2 n=6: implication holds; perms with no cut satisfying the criterion: 0
cor 5.2 n=7: implication holds; perms with no cut satisfying the criterion: 0
cor 5.2 n=8: implication holds; perms with no cut satisfying the criterion: 16, e.g. [0, 3, 6, 1, 4, 7, 2, 5]
verdict: PASS (28.0 s)
```
