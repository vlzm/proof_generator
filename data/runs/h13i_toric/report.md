# H13-I attempt (session 9): report

Version h13i_toric-1.0, core oracle-1.0, elapsed 91.9 s.

## Single-axis (A = 0, vary B only) -- REFUTED as sufficient

- n=4: max=2, bound=2, fails on 0 of n! perms, example (0, 2, 3, 1)
- n=5: max=4, bound=4, fails on 0 of n! perms, example (0, 2, 4, 3, 1)
- n=6: max=6, bound=6, fails on 0 of n! perms, example (0, 2, 4, 5, 3, 1)
- n=7: max=10, bound=9, fails on 7 of n! perms, example (0, 5, 4, 3, 2, 1, 6)
- n=8: max=13, bound=12, fails on 16 of n! perms, example (0, 6, 5, 3, 4, 2, 1, 7)

## Diagonal + antidiagonal (2n candidates) -- REFUTED as sufficient

- n=4: max=3, bound=2, fails on 4 of n! perms, example (0, 3, 2, 1)
- n=5: max=5, bound=4, fails on 6 of n! perms, example (0, 3, 4, 2, 1)
- n=6: max=8, bound=6, fails on 26 of n! perms, example (0, 4, 5, 3, 1, 2)
- n=7: max=10, bound=9, fails on 30 of n! perms, example (0, 5, 4, 6, 3, 1, 2)
- n=8: max=16, bound=12, fails on 427 of n! perms, example (3, 6, 5, 0, 7, 2, 1, 4)

## Exact block window of size floor(n/2) -- does NOT always exist

- n=4: k=2, perms without such a window: 0, example None
- n=5: k=2, perms without such a window: 10, example (0, 2, 4, 1, 3)
- n=6: k=3, perms without such a window: 204, example (0, 1, 3, 5, 2, 4)
- n=7: k=3, perms without such a window: 1694, example (0, 1, 3, 4, 6, 2, 5)
- n=8: k=4, perms without such a window: 24464, example (0, 1, 2, 4, 5, 7, 3, 6)

## Single-element induction -- REFUTED

- n=4: increment=1, fails on 0 perms, worst gap None, example None
- n=5: increment=2, fails on 0 perms, worst gap None, example None
- n=6: increment=2, fails on 0 perms, worst gap None, example None
- n=7: increment=3, fails on 0 perms, worst gap None, example None
- n=8: increment=3, fails on 16 perms, worst gap 2, example (0, 5, 2, 7, 4, 1, 6, 3)

## Pair-element induction (best pair removed) -- REFUTED

- n=4: increment=2, fails on 0 perms, worst gap None, example None
- n=5: increment=3, fails on 0 perms, worst gap None, example None
- n=6: increment=4, fails on 0 perms, worst gap None, example None
- n=7: increment=5, fails on 0 perms, worst gap None, example None
- n=8: increment=6, fails on 8 perms, worst gap 1, example (0, 5, 2, 7, 4, 1, 6, 3)

## Reflections: I(pi_h) = floor((n-1)^2/4) exactly, all h -- PROVED (all n) + VERIFIED

- n=4: bound=2, all h match: True
- n=5: bound=4, all h match: True
- n=6: bound=6, all h match: True
- n=7: bound=9, all h match: True
- n=8: bound=12, all h match: True
- n=9: bound=16, all h match: True
- n=10: bound=20, all h match: True
- n=11: bound=25, all h match: True
- n=12: bound=30, all h match: True
- n=13: bound=36, all h match: True
- n=14: bound=42, all h match: True
- n=15: bound=49, all h match: True
- n=16: bound=56, all h match: True
- n=17: bound=64, all h match: True
- n=18: bound=72, all h match: True
- n=19: bound=81, all h match: True
- n=20: bound=90, all h match: True
- n=21: bound=100, all h match: True
- n=22: bound=110, all h match: True
- n=23: bound=121, all h match: True
- n=24: bound=132, all h match: True
- n=25: bound=144, all h match: True
- n=26: bound=156, all h match: True
- n=27: bound=169, all h match: True
- n=28: bound=182, all h match: True
- n=29: bound=196, all h match: True
- n=30: bound=210, all h match: True
- n=31: bound=225, all h match: True
- n=32: bound=240, all h match: True
- n=33: bound=256, all h match: True
- n=34: bound=272, all h match: True
- n=35: bound=289, all h match: True
- n=36: bound=306, all h match: True
- n=37: bound=324, all h match: True
- n=38: bound=342, all h match: True
- n=39: bound=361, all h match: True
- n=40: bound=380, all h match: True

## Elementary proof for reflections

pi_h(i) = (h - i) mod n. Cut A = 0, B = h - m (m = 0..n-1 free): w_j = (m - j) mod n = m, m-1, ..., 0, n-1, n-2, ..., m+1. Two descending blocks of sizes m+1 and n-1-m; every cross-block pair is concordant (block 1 values 0..m all smaller than block 2 values m+1..n-1, and block 1 comes first), so inv(w) = C(m+1,2) + C(n-1-m,2) exactly. Minimising the convex function of m over m=0..n-1 gives exactly floor((n-1)^2/4) (matches C12's B_n = floor(n^2/4) + floor((n-1)^2/4) split). No case split on n mod 4 or on h is needed; this holds for every n >= 1.
