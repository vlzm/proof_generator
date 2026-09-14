# H13-I structured/random sampling -- report

Version h13i_sample-1.0, core oracle-1.0, seed 42, 21.4 s total.
SAMPLED coverage only (not exhaustive); exhaustive range is h13i_toric.c (n <= 12).

| n | target | worst observed |
|---|---|---|
| 11 | 25 | 25 |
| 12 | 30 | 30 |
| 15 | 49 | 49 |
| 20 | 90 | 90 |
| 25 | 144 | 144 |
| 30 | 210 | 210 |
| 40 | 380 | 380 |
| 50 | 600 | 600 |
| 60 | 870 | 870 |
| 80 | 1560 | 1560 |
| 100 | 2450 | 2450 |
| 150 | 5550 | 5550 |
| 200 | 9900 | 9900 |

no counterexample found; reflections achieve I(pi) == target exactly at every n tested (consistent with C33's finite data: reflections are the extremal family)
