# toric_cuts — двойные разрезы: тождества и запрет усреднений

Команда: `python3 experiments/toric_cuts.py --nmax 8 --exhaustive-max 7 --sample 4000 --seed 20260915`. Версия: toric_cuts-1.0.

```text
== toric_cuts-1.0 args={'nmax': 8, 'exhaustive_max': 7, 'sample': 4000, 'seed': 20260915}
n=4 exhaustive (24 perms): max I = 2 (floor((n-1)^2/4) = 2) at (0, 3, 2, 1); identities pair=True weight=True triples=True; 0.0 s
  family uniform  : max_pi (sum over the best family) = 14 vs n*floor((n-1)^2/4) = 8 -> FAILS the bound; worst pi = (0, 3, 2, 1)
  family diag     : max_pi (sum over the best family) = 12 vs n*floor((n-1)^2/4) = 8 -> FAILS the bound; worst pi = (0, 3, 2, 1)
  family antidiag : max_pi (sum over the best family) = 8 vs n*floor((n-1)^2/4) = 8 -> within bound; worst pi = (0, 1, 2, 3)
  family slope    : max_pi (sum over the best family) = 8 vs n*floor((n-1)^2/4) = 8 -> within bound; worst pi = (0, 1, 3, 2)
  family anchored : max_pi (sum over the best family) = 8 vs n*floor((n-1)^2/4) = 8 -> within bound; worst pi = (0, 3, 2, 1)
n=5 exhaustive (120 perms): max I = 4 (floor((n-1)^2/4) = 4) at (0, 4, 3, 2, 1); identities pair=True weight=True triples=True; 0.0 s
  family uniform  : max_pi (sum over the best family) = 30 vs n*floor((n-1)^2/4) = 20 -> FAILS the bound; worst pi = (0, 4, 3, 2, 1)
  family diag     : max_pi (sum over the best family) = 30 vs n*floor((n-1)^2/4) = 20 -> FAILS the bound; worst pi = (0, 4, 3, 2, 1)
  family antidiag : max_pi (sum over the best family) = 25 vs n*floor((n-1)^2/4) = 20 -> FAILS the bound; worst pi = (0, 2, 4, 1, 3)
  family slope    : max_pi (sum over the best family) = 22 vs n*floor((n-1)^2/4) = 20 -> FAILS the bound; worst pi = (0, 2, 1, 4, 3)
  family anchored : max_pi (sum over the best family) = 21 vs n*floor((n-1)^2/4) = 20 -> FAILS the bound; worst pi = (0, 1, 4, 3, 2)
n=6 exhaustive (720 perms): max I = 6 (floor((n-1)^2/4) = 6) at (0, 5, 4, 3, 2, 1); identities pair=True weight=True triples=True; 0.1 s
  family uniform  : max_pi (sum over the best family) = 55 vs n*floor((n-1)^2/4) = 36 -> FAILS the bound; worst pi = (0, 5, 4, 3, 2, 1)
  family diag     : max_pi (sum over the best family) = 52 vs n*floor((n-1)^2/4) = 36 -> FAILS the bound; worst pi = (0, 5, 4, 3, 2, 1)
  family antidiag : max_pi (sum over the best family) = 44 vs n*floor((n-1)^2/4) = 36 -> FAILS the bound; worst pi = (0, 3, 2, 5, 4, 1)
  family slope    : max_pi (sum over the best family) = 38 vs n*floor((n-1)^2/4) = 36 -> FAILS the bound; worst pi = (0, 2, 5, 1, 4, 3)
  family anchored : max_pi (sum over the best family) = 36 vs n*floor((n-1)^2/4) = 36 -> within bound; worst pi = (0, 1, 4, 5, 3, 2)
n=7 exhaustive (5040 perms): max I = 9 (floor((n-1)^2/4) = 9) at (0, 6, 5, 4, 3, 2, 1); identities pair=True weight=True triples=True; 0.8 s
  family uniform  : max_pi (sum over the best family) = 91 vs n*floor((n-1)^2/4) = 63 -> FAILS the bound; worst pi = (0, 6, 5, 4, 3, 2, 1)
  family diag     : max_pi (sum over the best family) = 91 vs n*floor((n-1)^2/4) = 63 -> FAILS the bound; worst pi = (0, 6, 5, 4, 3, 2, 1)
  family antidiag : max_pi (sum over the best family) = 77 vs n*floor((n-1)^2/4) = 63 -> FAILS the bound; worst pi = (0, 3, 6, 2, 5, 1, 4)
  family slope    : max_pi (sum over the best family) = 68 vs n*floor((n-1)^2/4) = 63 -> FAILS the bound; worst pi = (0, 2, 4, 1, 6, 5, 3)
  family anchored : max_pi (sum over the best family) = 67 vs n*floor((n-1)^2/4) = 63 -> FAILS the bound; worst pi = (0, 2, 1, 6, 5, 4, 3)
n=8 sampled (4040 perms: random + affine + reflections): max I = 12 (floor((n-1)^2/4) = 12) at (5, 4, 3, 2, 1, 0, 7, 6); identities pair=True weight=True triples=True; 0.8 s
  family uniform  : max_pi (sum over the best family) = 140 vs n*floor((n-1)^2/4) = 96 -> FAILS the bound; worst pi = (5, 4, 3, 2, 1, 0, 7, 6)
  family diag     : max_pi (sum over the best family) = 136 vs n*floor((n-1)^2/4) = 96 -> FAILS the bound; worst pi = (5, 4, 3, 2, 1, 0, 7, 6)
  family antidiag : max_pi (sum over the best family) = 112 vs n*floor((n-1)^2/4) = 96 -> FAILS the bound; worst pi = (0, 7, 3, 6, 2, 1, 5, 4)
  family slope    : max_pi (sum over the best family) = 104 vs n*floor((n-1)^2/4) = 96 -> FAILS the bound; worst pi = (3, 6, 2, 1, 7, 5, 4, 0)
  family anchored : max_pi (sum over the best family) = 100 vs n*floor((n-1)^2/4) = 96 -> FAILS the bound; worst pi = (5, 7, 6, 4, 3, 0, 2, 1)

affine no-go (inv depends only on r = a*p + b - s):
  n=4: confirmed (8 affine inputs; every cut family meeting each residue r equally averages to mean_r inv)
  n=5: confirmed (20 affine inputs; every cut family meeting each residue r equally averages to mean_r inv)
  n=6: confirmed (12 affine inputs; every cut family meeting each residue r equally averages to mean_r inv)
  n=7: confirmed (42 affine inputs; every cut family meeting each residue r equally averages to mean_r inv)
  n=8: confirmed (32 affine inputs; every cut family meeting each residue r equally averages to mean_r inv)
  reflection n=4: mean inv over all cuts = 3.500 = C(n,2) - (n^2-1)/6 = 3.500; min = 2
  reflection n=5: mean inv over all cuts = 6.000 = C(n,2) - (n^2-1)/6 = 6.000; min = 4
  reflection n=6: mean inv over all cuts = 9.167 = C(n,2) - (n^2-1)/6 = 9.167; min = 6
  reflection n=7: mean inv over all cuts = 13.000 = C(n,2) - (n^2-1)/6 = 13.000; min = 9
  reflection n=8: mean inv over all cuts = 17.500 = C(n,2) - (n^2-1)/6 = 17.500; min = 12
verdict: PASS
```
