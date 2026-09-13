# Toric correlation identity — report

Version toric_corr-1.0; 13.3 s.

## Theorem A (identity), brute force over all cuts

- n = 2: 8 cuts, failures 0
- n = 3: 54 cuts, failures 0
- n = 4: 384 cuts, failures 0
- n = 5: 3000 cuts, failures 0
- n = 6: 25920 cuts, failures 0

## Theorem B (energy bound) and its equality set

- n = 2: bound 4, min E 4, below 0, equality on 2 perms, = reflections: True
- n = 3: bound 12, min E 12, below 0, equality on 3 perms, = reflections: True
- n = 4: bound 16, min E 16, below 0, equality on 4 perms, = reflections: True
- n = 5: bound 0, min E 0, below 0, equality on 5 perms, = reflections: True
- n = 6: bound -60, min E -60, below 0, equality on 6 perms, = reflections: True
- n = 7: bound -196, min E -196, below 0, equality on 7 perms, = reflections: True
- n = 8: bound -448, min E -448, below 0, equality on 8 perms, = reflections: True

## Corollary C: I(pi) <= floor((n-1)(2n-1)/6)

- n = 4: proved bound 3, max_pi I = 2, target floor((n-1)^2/4) = 2, violations 0
- n = 5: proved bound 6, max_pi I = 4, target floor((n-1)^2/4) = 4, violations 0
- n = 6: proved bound 9, max_pi I = 6, target floor((n-1)^2/4) = 6, violations 0
- n = 7: proved bound 13, max_pi I = 9, target floor((n-1)^2/4) = 9, violations 0
- n = 8: proved bound 17, max_pi I = 12, target floor((n-1)^2/4) = 12, violations 0

## Candidate sufficient conditions for H13-I (all fail)

- n = 4 (24 perms), thr = 14, reflection peak = 12: failures L=0, M=20, P=0, T=24, W=0, exact=0
    first M-failure: [0, 1, 3, 2]
    first T-failure: [0, 1, 2, 3]
- n = 5 (120 perms), thr = 20, reflection peak = 20: failures L=25, M=90, P=0, T=0, W=25, exact=0
    first L-failure: [0, 2, 1, 4, 3]
    first M-failure: [0, 1, 3, 4, 2]
    first W-failure: [0, 1, 4, 3, 2]
- n = 6 (720 perms), thr = 33, reflection peak = 38: failures L=72, M=540, P=0, T=720, W=0, exact=0
    first L-failure: [0, 2, 5, 1, 4, 3]
    first M-failure: [0, 1, 2, 5, 4, 3]
    first T-failure: [0, 1, 2, 3, 4, 5]
- n = 7 (5040 perms), thr = 42, reflection peak = 56: failures L=343, M=2583, P=0, T=0, W=441, exact=0
    first L-failure: [0, 2, 4, 1, 6, 5, 3]
    first M-failure: [0, 1, 2, 5, 6, 3, 4]
    first W-failure: [0, 1, 4, 6, 5, 3, 2]
- n = 8 (40320 perms), thr = 60, reflection peak = 88: failures L=1152, M=20456, P=16, T=40296, W=1088, exact=0
    first L-failure: [0, 2, 6, 5, 1, 4, 7, 3]
    first M-failure: [0, 1, 2, 5, 6, 7, 3, 4]
    first P-failure: [0, 3, 6, 1, 4, 7, 2, 5]
    first T-failure: [0, 1, 2, 3, 4, 5, 6, 7]
    first W-failure: [0, 1, 4, 7, 6, 5, 3, 2]

## Peak bound on random permutations

- n = 12: reflection peak 292, min random max S 320, fails: False
- n = 20: reflection peak 1340, min random max S 1212, fails: True
- n = 30: reflection peak 4510, min random max S 3710, fails: True
- n = 50: reflection peak 20850, min random max S 13618, fails: True
- n = 80: reflection peak 85360, min random max S 45052, fails: True
