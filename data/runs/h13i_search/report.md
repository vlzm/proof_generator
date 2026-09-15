# H13-I candidate search (session 9)

Version h13i_search-1.0. Exhaustive over all pi, 4 <= n <= 8.

## Candidate A: fixed diagonal k = (a-b) mod n

REFUTED for every n and every k: worst-case min_a f(a,a-k) exceeds floor((n-1)^2/4) for every choice of k.

| n | bound | worst any k | global I_max (full n^2 search) |
|---|---|---|---|
| 4 | 2 | 4 | 2 |
| 5 | 4 | 6 | 4 |
| 6 | 6 | 10 | 6 |
| 7 | 9 | 12 | 9 |
| 8 | 12 | 17 | 12 |

## Candidate B: adjacent-optimal-b chain (jump <= 1)

REFUTED: fails for the majority of pi already at n = 5.

| n | pi without short chain | total pi |
|---|---|---|
| 4 | 8 | 24 |
| 5 | 80 | 120 |
| 6 | 660 | 720 |
| 7 | 4725 | 5040 |
| 8 | 39824 | 40320 |

Example pi with no jump<=1 chain (n = smallest tested): [[0, 2, 1, 3], [0, 3, 1, 2], [1, 0, 2, 3], [1, 3, 2, 0], [2, 0, 3, 1]]
