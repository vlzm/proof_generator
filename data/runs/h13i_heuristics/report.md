# h13i_heuristics — три неудачные эвристики выбора (a, b) для H13-I

Команда: `python3 experiments/h13i_heuristics.py --nmax 8`. Версия: h13i_heuristics-1.0.

```text
== h13i_heuristics-1.0 args={'nmax': 8}
n=4 H1_fixed_a0: worst over all pi = 2 (bound 2) within bound; witness (0, 2, 3, 1)
n=4 H2_greedy_two_step: worst over all pi = 2 (bound 2) within bound; witness (0, 3, 2, 1)
n=4 H3_canonical_b: worst over all pi = 3 (bound 2) EXCEEDS; witness (0, 3, 2, 1)
n=4: 0.00 s
n=5 H1_fixed_a0: worst over all pi = 4 (bound 4) within bound; witness (0, 2, 4, 3, 1)
n=5 H2_greedy_two_step: worst over all pi = 4 (bound 4) within bound; witness (0, 4, 3, 2, 1)
n=5 H3_canonical_b: worst over all pi = 6 (bound 4) EXCEEDS; witness (0, 4, 3, 2, 1)
n=5: 0.00 s
n=6 H1_fixed_a0: worst over all pi = 6 (bound 6) within bound; witness (0, 2, 4, 5, 3, 1)
n=6 H2_greedy_two_step: worst over all pi = 6 (bound 6) within bound; witness (0, 4, 3, 2, 1, 5)
n=6 H3_canonical_b: worst over all pi = 10 (bound 6) EXCEEDS; witness (0, 5, 4, 3, 2, 1)
n=6: 0.04 s
n=7 H1_fixed_a0: worst over all pi = 10 (bound 9) EXCEEDS; witness (0, 5, 4, 3, 2, 1, 6)
n=7 H2_greedy_two_step: worst over all pi = 10 (bound 9) EXCEEDS; witness (0, 5, 4, 3, 2, 1, 6)
n=7 H3_canonical_b: worst over all pi = 15 (bound 9) EXCEEDS; witness (0, 6, 5, 4, 3, 2, 1)
n=7: 0.31 s
n=8 H1_fixed_a0: worst over all pi = 13 (bound 12) EXCEEDS; witness (0, 6, 5, 3, 4, 2, 1, 7)
n=8 H2_greedy_two_step: worst over all pi = 12 (bound 12) within bound; witness (0, 5, 4, 3, 2, 7, 1, 6)
n=8 H3_canonical_b: worst over all pi = 21 (bound 12) EXCEEDS; witness (0, 7, 6, 5, 4, 3, 2, 1)
n=8: 3.39 s
```
