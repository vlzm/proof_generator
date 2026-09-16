# H13-I reduction probe (session 9)

Version h13i_reduction_probe-1.0, core oracle-1.0.

| n | bound | max I | corner fails | pair-removal budget | max gap | pair-removal fails |
|---|---|---|---|---|---|---|
| 4 | 2 | 2 | 4 | 2 | 2 | 0 |
| 5 | 4 | 4 | 5 | 3 | 3 | 0 |
| 6 | 6 | 6 | 6 | 4 | 4 | 0 |
| 7 | 9 | 9 | 56 | 5 | 5 | 0 |
| 8 | 12 | 12 | 552 | 6 | 7 | 8 |

Conclusion: corner cuts (n candidates through a data point) fail already at n=4; two-point removal induction (budget n-2) holds exactly for 4<=n<=7 (equality on reflections) and fails at n=8 on exactly the 8 rotations of pi(i) = 5i mod 8 (gap 7 > budget 6). Both approaches are ruled out as proof strategies for H13-I.
