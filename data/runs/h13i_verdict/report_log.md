# toric_cut_algebra — лог прогона

Команда: `python3 experiments/toric_cut_algebra.py --nmax 8 --free-nmax 8`. Версии: toric_cut-1.0, oracle-1.0.

```text
== toric_cut-1.0 oracle-1.0 args={'nmax': 8, 'free_nmax': 8}
part 1 n=4: 24 perms, (q,c) and (m,k) families of lines coincide -> ok
part 1 n=5: 120 perms, (q,c) and (m,k) families of lines coincide -> ok
part 1 n=6: 720 perms, (q,c) and (m,k) families of lines coincide -> ok
part 2 n=4: 24 perms x 16 cuts, identity L1 holds -> ok (0.0 s)
part 2 n=5: 120 perms x 25 cuts, identity L1 holds -> ok (0.0 s)
part 2 n=6: 720 perms x 36 cuts, identity L1 holds -> ok (0.1 s)
part 2 n=7: 5040 perms x 49 cuts, identity L1 holds -> ok (1.6 s)
part 2 n=8: 200 random (pi, m, k), identity L1 holds -> ok
part 2 n=12: 200 random (pi, m, k), identity L1 holds -> ok
part 2 n=20: 200 random (pi, m, k), identity L1 holds -> ok
part 2 n=33: 200 random (pi, m, k), identity L1 holds -> ok
part 2 n=50: 200 random (pi, m, k), identity L1 holds -> ok
part 3 n=4: pair formula (n-d_p)(n-d_v)+d_p d_v and Sigma(pi) -> ok
part 3 n=5: pair formula (n-d_p)(n-d_v)+d_p d_v and Sigma(pi) -> ok
part 3 n=6: pair formula (n-d_p)(n-d_v)+d_p d_v and Sigma(pi) -> ok
part 4 n=4: target floor(n^2/4)=4; mean(id)=3.500 (= (n-1)(2n-1)/6), mean(sigma_n)=2.500 (= (n^2-1)/6); averaging settles 0/24 perms
part 4 n=5: target floor(n^2/4)=6; mean(id)=6.000 (= (n-1)(2n-1)/6), mean(sigma_n)=4.000 (= (n^2-1)/6); averaging settles 5/120 perms
part 4 n=6: target floor(n^2/4)=9; mean(id)=9.167 (= (n-1)(2n-1)/6), mean(sigma_n)=5.833 (= (n^2-1)/6); averaging settles 6/720 perms
part 4 n=7: target floor(n^2/4)=12; mean(id)=13.000 (= (n-1)(2n-1)/6), mean(sigma_n)=8.000 (= (n^2-1)/6); averaging settles 56/5040 perms
part 4 n=8: target floor(n^2/4)=16; mean(id)=17.500 (= (n-1)(2n-1)/6), mean(sigma_n)=10.500 (= (n^2-1)/6); averaging settles 488/40320 perms
part 5 n=4: max_pi min_free = 2 (target floor((n-1)^2/4) = 2), max_pi I(pi) = 2; I - min_free histogram {0: 24}; of the 0 perms with a gap 0 are affine -> ok
part 5 n=5: max_pi min_free = 4 (target floor((n-1)^2/4) = 4), max_pi I(pi) = 4; I - min_free histogram {0: 120}; of the 0 perms with a gap 0 are affine -> ok
part 5 n=6: max_pi min_free = 6 (target floor((n-1)^2/4) = 6), max_pi I(pi) = 6; I - min_free histogram {0: 720}; of the 0 perms with a gap 0 are affine -> ok
part 5 n=7: max_pi min_free = 9 (target floor((n-1)^2/4) = 9), max_pi I(pi) = 9; I - min_free histogram {0: 5040}; of the 0 perms with a gap 0 are affine -> ok
part 5 n=8: max_pi min_free = 12 (target floor((n-1)^2/4) = 12), max_pi I(pi) = 12; I - min_free histogram {0: 40112, 1: 200, 3: 8}; of the 208 perms with a gap 16 are affine -> ok
part 5 n=12: local search over free subsets on 43 inputs, max = 30 <= 30 -> ok
part 5 n=20: local search over free subsets on 63 inputs, max = 90 <= 90 -> ok
part 5 n=33: local search over free subsets on 112 inputs, max = 256 <= 256 -> ok
part 5 n=50: local search over free subsets on 129 inputs, max = 600 <= 600 -> ok
part 5 n=100: local search over free subsets on 239 inputs, max = 2450 <= 2450 -> ok
part 5 n=101: local search over free subsets on 420 inputs, max = 2500 <= 2500 -> ok
part 6 n=4: cuts failing criterion H14 for 0 perms
part 6 n=5: cuts failing criterion H14 for 0 perms
part 6 n=6: cuts failing criterion H14 for 0 perms
part 6 n=7: cuts failing criterion H14 for 0 perms
part 6 n=8: cuts failing criterion H14 for 16 perms (16 affine); e.g. [0, 3, 6, 1, 4, 7, 2, 5]
part 7 n=4: max_pi I(pi) = 2 (target 2), attained on 4 perms -> ok
part 7 n=5: max_pi I(pi) = 4 (target 4), attained on 5 perms -> ok
part 7 n=6: max_pi I(pi) = 6 (target 6), attained on 6 perms -> ok
part 7 n=7: max_pi I(pi) = 9 (target 9), attained on 7 perms -> ok
part 7 n=8: max_pi I(pi) = 12 (target 12), attained on 8 perms -> ok
part 7 n=20: 63 structured/random inputs, max I = 90 <= 90 -> ok
part 7 n=50: 75 structured/random inputs, max I = 600 <= 600 -> ok
part 7 n=100: 64 structured/random inputs, max I = 2450 <= 2450 -> ok
part 7 n=101: 60 structured/random inputs, max I = 2500 <= 2500 -> ok
verdict: PASS (82.2 s)
```
