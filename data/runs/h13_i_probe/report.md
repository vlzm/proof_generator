# H13-I probe (h13_i_probe-1.0)

Seed 0, 36 s.

## Multiplicative family pi(i) = k*i mod n (all k coprime to n)

| n | floor((n-1)^2/4) | max I over k | argmax k | ratio |
|---|---|---|---|---|
| 5 | 4 | 4 | 4 | 1.0 |
| 6 | 6 | 6 | 5 | 1.0 |
| 7 | 9 | 9 | 6 | 1.0 |
| 8 | 12 | 12 | 7 | 1.0 |
| 9 | 16 | 16 | 8 | 1.0 |
| 10 | 20 | 20 | 9 | 1.0 |
| 11 | 25 | 25 | 10 | 1.0 |
| 12 | 30 | 30 | 11 | 1.0 |
| 13 | 36 | 36 | 12 | 1.0 |
| 14 | 42 | 42 | 13 | 1.0 |
| 15 | 49 | 49 | 14 | 1.0 |
| 16 | 56 | 56 | 15 | 1.0 |
| 17 | 64 | 64 | 16 | 1.0 |
| 18 | 72 | 72 | 17 | 1.0 |
| 19 | 81 | 81 | 18 | 1.0 |
| 20 | 90 | 90 | 19 | 1.0 |
| 21 | 100 | 100 | 20 | 1.0 |
| 22 | 110 | 110 | 21 | 1.0 |
| 23 | 121 | 121 | 22 | 1.0 |
| 24 | 132 | 132 | 23 | 1.0 |
| 25 | 144 | 144 | 24 | 1.0 |
| 26 | 156 | 156 | 25 | 1.0 |
| 27 | 169 | 169 | 26 | 1.0 |
| 28 | 182 | 182 | 27 | 1.0 |
| 29 | 196 | 196 | 28 | 1.0 |
| 30 | 210 | 210 | 29 | 1.0 |
| 31 | 225 | 225 | 30 | 1.0 |

Worst ratio among non-reflection multipliers (k != n-1): {'n': 29, 'floor': 196, 'max_I_over_k': 196, 'argmax_k': 28, 'ratio': 1.0, 'best_nonreflection_k': 19, 'best_nonreflection_I': 191, 'best_nonreflection_ratio': 0.9745}

## Hill-climbing / reflection-neighborhood check

| n | floor | I(reflection) | best after 1 transposition | hill-climb best |
|---|---|---|---|---|
| 12 | 30 | 30 | 30 | 25 |
| 15 | 49 | 49 | 49 | 42 |
| 18 | 72 | 72 | 72 | 71 |
| 20 | 90 | 90 | 90 | 79 |
| 24 | 132 | 132 | 132 | 115 |

## Sliding-window sub-lemma

FALSE in general: n=5, pi(i) = 2*i mod 5, no window of size 2 maps onto a contiguous value-window (checked all n windows).

## Conclusion

no violation of H13-I (I(pi) <= floor((n-1)^2/4)) found; reflection sigma_n is the unique maximizer within tested families and a strict local max under single transpositions up to n=24.
