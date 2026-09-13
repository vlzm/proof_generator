# check_C27_reflections — summary

- odd n 5..60: F_c=P and S_c=3(n-1)/2 at every (h,c) -- 35980 triples; 2F_1-S_1+H(1)=B_n exactly at every h -- 896 cases
- even n 4..60: closed-form F_c,S_c,H(c) match construction at c=0,1 for every h -- 1856 (h,c) pairs; excess pairs by class: n%4=0,h%2=0:(3, -1), n%4=0,h%2=1:(0, 2), n%4=2,h%2=0:(1, 1), n%4=2,h%2=1:(2, 0)
- route=n1: all reflections 4<=n<=40, all shifts, words executed and verified to sort (22126 words); worst min_c(len)-B_n = 1 at (6, 0)
- route=carrier: all reflections 4<=n<=40, all shifts, words executed and verified to sort (22126 words); worst min_c(len)-B_n = 1 at (6, 2)

Total time: 45.4 s.
