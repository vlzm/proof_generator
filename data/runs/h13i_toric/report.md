# H13-I toric reformulation: identity check + three refuted families

Version h13i_toric-1.0, elapsed 0.4s.

## Pairwise inversion-count identity

`sum_{a,b} inv(w^{a,b}) = sum_{i<i'} [n(d+e) - 2de]` verified by direct
computation (random permutations, `--formula-trials` each) for all n in
range checked (3..9): all OK.

## Weakened cut families (exhaustive over all pi, per n, up to the smallest
failing n or `--amax`)

### single_a (b=0 fixed)

| n | target floor((n-1)^2/4) | max over pi | sufficient? | worst pi |
|---|---|---|---|---|
| 3 | 1 | 1 | True | None |
| 4 | 2 | 2 | True | None |
| 5 | 4 | 4 | True | None |
| 6 | 6 | 6 | True | None |
| 7 | 9 | 10 | False | [0, 5, 4, 3, 2, 1, 6] |

### a_then_avg (best a, pigeonhole on b)

| n | target floor((n-1)^2/4) | max over pi | sufficient? | worst pi |
|---|---|---|---|---|
| 3 | 1 | 1 | True | None |
| 4 | 2 | 3 | False | [0, 1, 3, 2] |

### vertex_cut (best data point as origin)

| n | target floor((n-1)^2/4) | max over pi | sufficient? | worst pi |
|---|---|---|---|---|
| 3 | 1 | 1 | True | None |
| 4 | 2 | 3 | False | [0, 3, 2, 1] |

