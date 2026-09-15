# h3_embed_probe report (h3_embed_probe-1.0, oracle-1.0)

Naive embedding: append one element at the end of the array, replay the same
move string that sorts sigma_n0 verbatim on the (n0+1)-array, measure the
distance from the result to sigma_n.

| n0 | n | len(word) | N_rot(word) | patch cost | budget n-1 | within budget? | total | B_n |
|---|---|---|---|---|---|---|---|---|
| 4 | 5 | 6 | 4 | 4 | 4 | True | 10 | 10 |
| 5 | 6 | 10 | 6 | 7 | 5 | False | 17 | 15 |
| 6 | 7 | 15 | 9 | 7 | 6 | False | 22 | 21 |
| 7 | 8 | 21 | 12 | 7 | 7 | True | 28 | 28 |
| 8 | 9 | 28 | 16 | 10 | 8 | False | 38 | 36 |
| 9 | 10 | 36 | 20 | 13 | 9 | False | 49 | 45 |

4/6 probes exceed the n-1 budget -> naive replay is not the H3 mechanism.
