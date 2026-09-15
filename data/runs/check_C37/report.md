# check_C37 — торическая переформулировка H13-I

Команда: `python3 checks/check_C37.py --amax 8`. Версии: check_C37-1.0, oracle-1.0.

```text
== check_C37-1.0 oracle-1.0 args={'amax': 8}
part A n=4: 4608 (pi, cut, ordered pair) checks, 0.0 s
part A n=5: 60000 (pi, cut, ordered pair) checks, 0.0 s
part A n=6: 777600 (pi, cut, ordered pair) checks, 0.4 s
part B n=4: 384 cuts, gradient ok, 0.0 s
part B n=5: 3000 cuts, gradient ok, 0.0 s
part B n=6: 25920 cuts, gradient ok, 0.1 s
part B n=7: 246960 cuts, gradient ok, 0.9 s
part C n=4: 6144 (cut, shift) pairs, block shift ok, 0.0 s
part C n=5: 75000 (cut, shift) pairs, block shift ok, 0.2 s
part D n=4: 14 123-avoiding lines, inv >= 2, 1 with equality = reflection lines, 0.0 s
part D n=5: 42 123-avoiding lines, inv >= 4, 2 with equality = reflection lines, 0.0 s
part D n=6: 132 123-avoiding lines, inv >= 6, 1 with equality = reflection lines, 0.0 s
part D n=7: 429 123-avoiding lines, inv >= 9, 2 with equality = reflection lines, 0.0 s
part D n=8: 1430 123-avoiding lines, inv >= 12, 1 with equality = reflection lines, 0.0 s
part E n=4: 384 cuts, ascending triples are positive triples, T+=0 on 4 reflections, 0.0 s
part E n=5: 3000 cuts, ascending triples are positive triples, T+=0 on 5 reflections, 0.0 s
part E n=6: 25920 cuts, ascending triples are positive triples, T+=0 on 6 reflections, 0.1 s
part F n=4: max_pi I = 2 = floor((n-1)^2/4), attained on 4 reflections, 0.0 s
part F n=5: max_pi I = 4 = floor((n-1)^2/4), attained on 5 reflections, 0.0 s
part F n=6: max_pi I = 6 = floor((n-1)^2/4), attained on 6 reflections, 0.0 s
part F n=7: max_pi I = 9 = floor((n-1)^2/4), attained on 7 reflections, 0.1 s
part F n=8: max_pi I = 12 = floor((n-1)^2/4), attained on 8 reflections, 0.6 s
total 2.4 s
verdict: PASS
```
