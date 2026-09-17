# check_C37 — среднее по разрезам торического класса (леммы 1–7, теорема 3)

Команда: `python3 checks/check_C37.py --nmax1 6 --nmax2 8 --nmax4 24 --nmax5 7`. Версия: check_C37-1.0.
Доказательство: `docs/proofs/C37_toric_average.md`.

```text
== check_C37-1.0 args={'nmax1': 6, 'nmax2': 8, 'nmax4': 24, 'nmax5': 7}
part 1 (Lemmas 1, 7): n=4, 24 permutations, definition = pair count = incremental rule, 0.0 s
part 1 (Lemmas 1, 7): n=5, 120 permutations, definition = pair count = incremental rule, 0.0 s
part 1 (Lemmas 1, 7): n=6, 720 permutations, definition = pair count = incremental rule, 0.1 s
part 2-3 (Lemma 2, Theorem 3, Corollary 4): n=4, 24 permutations, identity exact, I <= 3 (target 2); Corollary 4 covers 0/24; Lemma 5: 4 permutations with Delta = 0 = the n reflections, 0.0 s
part 2-3 (Lemma 2, Theorem 3, Corollary 4): n=5, 120 permutations, identity exact, I <= 6 (target 4); Corollary 4 covers 0/120; Lemma 5: 5 permutations with Delta = 0 = the n reflections, 0.0 s
part 2-3 (Lemma 2, Theorem 3, Corollary 4): n=6, 720 permutations, identity exact, I <= 9 (target 6); Corollary 4 covers 6/720; Lemma 5: 6 permutations with Delta = 0 = the n reflections, 0.0 s
part 2-3 (Lemma 2, Theorem 3, Corollary 4): n=7, 5040 permutations, identity exact, I <= 13 (target 9); Corollary 4 covers 56/5040; Lemma 5: 7 permutations with Delta = 0 = the n reflections, 0.1 s
part 2-3 (Lemma 2, Theorem 3, Corollary 4): n=8, 40320 permutations, identity exact, I <= 17 (target 12); Corollary 4 covers 488/40320; Lemma 5: 8 permutations with Delta = 0 = the n reflections, 0.9 s
part 4 (Lemma 6): n=4, all h: inv = C(K+1,2)+C(n-K-1,2) at all 16 cuts, I(pi_h) = 2 = floor((n-1)^2/4) at 4 cuts
part 4 (Lemma 6): n=5, all h: inv = C(K+1,2)+C(n-K-1,2) at all 25 cuts, I(pi_h) = 4 = floor((n-1)^2/4) at 10 cuts
part 4 (Lemma 6): n=6, all h: inv = C(K+1,2)+C(n-K-1,2) at all 36 cuts, I(pi_h) = 6 = floor((n-1)^2/4) at 6 cuts
part 4 (Lemma 6): n=7, all h: inv = C(K+1,2)+C(n-K-1,2) at all 49 cuts, I(pi_h) = 9 = floor((n-1)^2/4) at 14 cuts
part 4 (Lemma 6): n=8, all h: inv = C(K+1,2)+C(n-K-1,2) at all 64 cuts, I(pi_h) = 12 = floor((n-1)^2/4) at 8 cuts
part 4 (Lemma 6): n=9, all h: inv = C(K+1,2)+C(n-K-1,2) at all 81 cuts, I(pi_h) = 16 = floor((n-1)^2/4) at 18 cuts
part 4 (Lemma 6): n=10, all h: inv = C(K+1,2)+C(n-K-1,2) at all 100 cuts, I(pi_h) = 20 = floor((n-1)^2/4) at 10 cuts
part 4 (Lemma 6): n=11, all h: inv = C(K+1,2)+C(n-K-1,2) at all 121 cuts, I(pi_h) = 25 = floor((n-1)^2/4) at 22 cuts
part 4 (Lemma 6): n=12, all h: inv = C(K+1,2)+C(n-K-1,2) at all 144 cuts, I(pi_h) = 30 = floor((n-1)^2/4) at 12 cuts
part 4 (Lemma 6): n=13, all h: inv = C(K+1,2)+C(n-K-1,2) at all 169 cuts, I(pi_h) = 36 = floor((n-1)^2/4) at 26 cuts
part 4 (Lemma 6): n=14, all h: inv = C(K+1,2)+C(n-K-1,2) at all 196 cuts, I(pi_h) = 42 = floor((n-1)^2/4) at 14 cuts
part 4 (Lemma 6): n=15, all h: inv = C(K+1,2)+C(n-K-1,2) at all 225 cuts, I(pi_h) = 49 = floor((n-1)^2/4) at 30 cuts
part 4 (Lemma 6): n=16, all h: inv = C(K+1,2)+C(n-K-1,2) at all 256 cuts, I(pi_h) = 56 = floor((n-1)^2/4) at 16 cuts
part 4 (Lemma 6): n=17, all h: inv = C(K+1,2)+C(n-K-1,2) at all 289 cuts, I(pi_h) = 64 = floor((n-1)^2/4) at 34 cuts
part 4 (Lemma 6): n=18, all h: inv = C(K+1,2)+C(n-K-1,2) at all 324 cuts, I(pi_h) = 72 = floor((n-1)^2/4) at 18 cuts
part 4 (Lemma 6): n=19, all h: inv = C(K+1,2)+C(n-K-1,2) at all 361 cuts, I(pi_h) = 81 = floor((n-1)^2/4) at 38 cuts
part 4 (Lemma 6): n=20, all h: inv = C(K+1,2)+C(n-K-1,2) at all 400 cuts, I(pi_h) = 90 = floor((n-1)^2/4) at 20 cuts
part 4 (Lemma 6): n=21, all h: inv = C(K+1,2)+C(n-K-1,2) at all 441 cuts, I(pi_h) = 100 = floor((n-1)^2/4) at 42 cuts
part 4 (Lemma 6): n=22, all h: inv = C(K+1,2)+C(n-K-1,2) at all 484 cuts, I(pi_h) = 110 = floor((n-1)^2/4) at 22 cuts
part 4 (Lemma 6): n=23, all h: inv = C(K+1,2)+C(n-K-1,2) at all 529 cuts, I(pi_h) = 121 = floor((n-1)^2/4) at 46 cuts
part 4 (Lemma 6): n=24, all h: inv = C(K+1,2)+C(n-K-1,2) at all 576 cuts, I(pi_h) = 132 = floor((n-1)^2/4) at 24 cuts
part 5 (H13-I finite range): n=4, 6 class representatives, max I = 2 = floor((n-1)^2/4) = 2, attained by 1 toric class, independent recomputation 2; min number of cuts with inv <= target = 4 -> ok
part 5 (H13-I finite range): n=5, 24 class representatives, max I = 4 = floor((n-1)^2/4) = 4, attained by 1 toric class, independent recomputation 4; min number of cuts with inv <= target = 7 -> ok
part 5 (H13-I finite range): n=6, 120 class representatives, max I = 6 = floor((n-1)^2/4) = 6, attained by 1 toric class, independent recomputation 6; min number of cuts with inv <= target = 6 -> ok
part 5 (H13-I finite range): n=7, 720 class representatives, max I = 9 = floor((n-1)^2/4) = 9, attained by 1 toric class, independent recomputation 9; min number of cuts with inv <= target = 10 -> ok
part 5 (H13-I finite range): n=8, 5040 class representatives, max I = 12 = floor((n-1)^2/4) = 12, attained by 1 toric class; min number of cuts with inv <= target = 8 -> ok
part 5 (H13-I finite range): n=9, 40320 class representatives, max I = 16 = floor((n-1)^2/4) = 16, attained by 1 toric class; min number of cuts with inv <= target = 14 -> ok
part 5 (H13-I finite range): n=10, 362880 class representatives, max I = 20 = floor((n-1)^2/4) = 20, attained by 1 toric class; min number of cuts with inv <= target = 10 -> ok
part 5 (H13-I finite range): n=11, 3628800 class representatives, max I = 25 = floor((n-1)^2/4) = 25, attained by 1 toric class; min number of cuts with inv <= target = 18 -> ok
part 5 (H13-I finite range): n=12, 39916800 class representatives, max I = 30 = floor((n-1)^2/4) = 30, attained by 1 toric class; min number of cuts with inv <= target = 12 -> ok
part 5 (H13-I finite range): n=13, 479001600 class representatives, max I = 36 = floor((n-1)^2/4) = 36, attained by 1 toric class; min number of cuts with inv <= target = 22 -> ok
verdict: PASS
```
