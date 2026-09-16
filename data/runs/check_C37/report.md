# check_C37 — проверка лемм доказательства C37/C38

Версия check_C37-1.0, seed 20260916, `--nbig 120`, 8 с. Результат: PASS.

Части 1-8 проверяют леммы доказательства (`docs/proofs/C37_toric_mean_bound.md`); часть 9 — только конечные
данные по открытой H13-I (не доказательство).

- part 1: lemma 2 and lemma 4 hold for n <= 2000
- part 2: lemma 1 verified on 924 permutations (all pi for n<=6, samples n=7..9)
- part 3-4: lemma 3 (exact mean) and lemma 5 (|A| bound, equality cases) hold
- part 5: C37 holds for all pi at 4<=n<=8 and on families up to n=120; max mean = (n-1)(2n-1)/6 exactly
- part 6: block criterion coverage n=4: 24/24, n=5: 110/120, n=6: 612/720, n=7: 3346/5040, n=8: 21616/40320; counterexample pi=2i mod 5 (block bound 6 > 4, I = 3)
- part 6b: averaging test (mean <= target, i.e. A >= n^3 for even n and n^2(n-1) for odd n) passes for n=4: 0/24, n=5: 5/120, n=6: 6/720, n=7: 56/5040, n=8: 488/40320
- part 7: lemma 6 / C38(b) hold for all reflections, 4 <= n <= 60 (optimal cuts: n for even n, 2n for odd n)
- part 8: the reformulation I <= floor((n-1)^2/4) <=> max D >= floor(n/2) verified for all pi, 4 <= n <= 7
- part 9 (evidence only): H13-I holds on all families and on the whole affine family for n up to 120; the maximum is always attained at a = n-1 (reflection)
