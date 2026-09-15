# toric_cut_peak: peak-rotation candidate for H13-I

Version: toric_cut_peak-1.0

| n | bound floor((n-1)^2/4) | worst (peak constr.) | fails | fails (global argmax only) | permutations | time (s) |
|---|---|---|---|---|---|---|
| 4 | 2 | 2 | 0 | 0 | 24 | 0.00 |
| 5 | 4 | 4 | 0 | 0 | 120 | 0.00 |
| 6 | 6 | 6 | 0 | 0 | 720 | 0.01 |
| 7 | 9 | 9 | 0 | 0 | 5040 | 0.08 |
| 8 | 12 | 12 | 0 | 0 | 40320 | 0.64 |
| 9 | 16 | 16 | 0 | 9 | 362880 | 8.07 |
| 10 | 20 | 20 | 0 | 20 | 3628800 | 81.67 |

peak-construction (all peaks, best of each) matches the bound exactly (no failures) for every pi at 4<=n<=10; the restricted global-argmax-only version fails at n in {9, 10} (see counterexamples); no general-n proof found this session

## Counterexamples to the global-argmax-only restriction

- n=9, pi=(0, 6, 5, 4, 3, 2, 8, 1, 7), bound=16: global-argmax peaks give [17, 17], but the full peak set gives [17, 15, 15, 17]
- n=9, pi=(1, 7, 0, 6, 5, 4, 3, 2, 8), bound=16: global-argmax peaks give [17, 17], but the full peak set gives [17, 17, 15, 15]
- n=9, pi=(2, 8, 1, 7, 0, 6, 5, 4, 3), bound=16: global-argmax peaks give [17, 17], but the full peak set gives [17, 17, 15, 15]
