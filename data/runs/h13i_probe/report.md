# H13-I probe (h13i_probe-1.0, core oracle-1.0)

n | bound | maxI | Laplacian | avg>bound | val-swap decr | pos-swap decr | time
---|---|---|---|---|---|---|---
4 | 2 | 2 | PASS | 24/24 | 32 | 32 | 0.01s
5 | 4 | 4 | PASS | 115/120 | 300 | 300 | 0.06s
6 | 6 | 6 | PASS | 714/720 | 1656 | 1656 | 0.76s
7 | 9 | 9 | PASS | 4984/5040 | 17640 | 17640 | 8.24s
8 | 12 | 12 | PASS | 39832/40320 | 134272 | 134272 | 111.35s

h13i_probe-1.0 core=oracle-1.0 nmax=8
n=4 bound=2 maxI=2 laplacian=PASS avg_exceeds=24/24 val_swap_decr=32 pos_swap_decr=32 (0.0s)
n=5 bound=4 maxI=4 laplacian=PASS avg_exceeds=115/120 val_swap_decr=300 pos_swap_decr=300 (0.1s)
n=6 bound=6 maxI=6 laplacian=PASS avg_exceeds=714/720 val_swap_decr=1656 pos_swap_decr=1656 (0.8s)
n=7 bound=9 maxI=9 laplacian=PASS avg_exceeds=4984/5040 val_swap_decr=17640 pos_swap_decr=17640 (8.2s)
n=8 bound=12 maxI=12 laplacian=PASS avg_exceeds=39832/40320 val_swap_decr=134272 pos_swap_decr=134272 (111.3s)
