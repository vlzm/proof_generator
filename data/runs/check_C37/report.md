# check_C37 report

Version check_C37-1.0. Overall: PASS

- part A: 19500 (n, pi, a) triples checked, ok=True
- part B: n=4 24 perms exhaustive, 0.0s, ok so far=True
- part B: n=5 120 perms exhaustive, 0.0s, ok so far=True
- part B: n=6 720 perms exhaustive, 0.1s, ok so far=True
- part B: n=7 5040 perms exhaustive, 0.7s, ok so far=True
- part B: n=8 40320 perms exhaustive, 8.6s, ok so far=True
- part B: 46224 permutations total, ok=True
- part C: n=4 max_I=2 == floor((n-1)^2/4)=2 (exhaustive, 24 perms)
- part C: n=5 max_I=4 == floor((n-1)^2/4)=4 (exhaustive, 120 perms)
- part C: n=6 max_I=6 == floor((n-1)^2/4)=6 (exhaustive, 720 perms)
- part C: n=7 max_I=9 == floor((n-1)^2/4)=9 (exhaustive, 5040 perms)
- part C: n=8 max_I=12 == floor((n-1)^2/4)=12 (exhaustive, 40320 perms)
- part C: n=9 max_I=16 == floor((n-1)^2/4)=16 (exhaustive, 362880 perms)
- part C: n=10 max_I=20 == floor((n-1)^2/4)=20 (exhaustive, 3628800 perms)
- part C: n=11 max_I=25 == floor((n-1)^2/4)=25 (exhaustive, 39916800 perms)
- part D: n=7 w=(0, 5, 4, 3, 2, 1, 6) inv(w)=10 min_b inv(rotate_b(w))=10 bound=9 (expected > bound, i.e. a-free strengthening refuted: True)
