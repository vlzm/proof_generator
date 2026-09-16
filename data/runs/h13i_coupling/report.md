
## Run 2026-09-16 11:05:54 (9.9 s)
args: {'nmax': 8, 'n': None, 'sample': None, 'seed': 1}

### n = 4 (exhaustive), bound floor((n-1)^2/4) = 2
- worst I(pi) (full n^2 search) = 2, example (0, 3, 2, 1)
- worst (fix a=0, vary b only) = 2, example (0, 2, 3, 1)
- worst (rule b=pi(a)) = 3, example (0, 3, 2, 1)
- worst (rule b=pi(a)-1) = 3, example (0, 1, 2, 3)
- worst (rule b=pi(a-1)+1) = 3, example (0, 3, 2, 1)

### n = 5 (exhaustive), bound floor((n-1)^2/4) = 4
- worst I(pi) (full n^2 search) = 4, example (0, 4, 3, 2, 1)
- worst (fix a=0, vary b only) = 4, example (0, 2, 4, 3, 1)
- worst (rule b=pi(a)) = 6, example (0, 4, 3, 2, 1)
- worst (rule b=pi(a)-1) = 5, example (0, 3, 1, 4, 2)
- worst (rule b=pi(a-1)+1) = 6, example (0, 4, 3, 2, 1)

### n = 6 (exhaustive), bound floor((n-1)^2/4) = 6
- worst I(pi) (full n^2 search) = 6, example (0, 5, 4, 3, 2, 1)
- worst (fix a=0, vary b only) = 6, example (0, 2, 4, 5, 3, 1)
- worst (rule b=pi(a)) = 10, example (0, 5, 4, 3, 2, 1)
- worst (rule b=pi(a)-1) = 7, example (0, 5, 4, 3, 2, 1)
- worst (rule b=pi(a-1)+1) = 10, example (0, 5, 4, 3, 2, 1)

### n = 7 (exhaustive), bound floor((n-1)^2/4) = 9
- worst I(pi) (full n^2 search) = 9, example (0, 6, 5, 4, 3, 2, 1)
- worst (fix a=0, vary b only) = 10, example (0, 5, 4, 3, 2, 1, 6)
- worst (rule b=pi(a)) = 15, example (0, 6, 5, 4, 3, 2, 1)
- worst (rule b=pi(a)-1) = 11, example (0, 5, 3, 1, 6, 4, 2)
- worst (rule b=pi(a-1)+1) = 15, example (0, 6, 5, 4, 3, 2, 1)

### n = 8 (exhaustive), bound floor((n-1)^2/4) = 12
- worst I(pi) (full n^2 search) = 12, example (0, 7, 6, 5, 4, 3, 2, 1)
- worst (fix a=0, vary b only) = 13, example (0, 6, 5, 3, 4, 2, 1, 7)
- worst (rule b=pi(a)) = 21, example (0, 7, 6, 5, 4, 3, 2, 1)
- worst (rule b=pi(a)-1) = 16, example (0, 7, 6, 5, 4, 3, 2, 1)
- worst (rule b=pi(a-1)+1) = 21, example (0, 7, 6, 5, 4, 3, 2, 1)

## Run 2026-09-16 11:07:44 (102.0 s)
args: {'nmax': None, 'n': 9, 'sample': None, 'seed': 1}

### n = 9 (exhaustive), bound floor((n-1)^2/4) = 16
- worst I(pi) (full n^2 search) = 16, example (0, 8, 7, 6, 5, 4, 3, 2, 1)
- worst (fix a=0, vary b only) = 17, example (0, 2, 7, 6, 5, 4, 3, 8, 1)
- worst (rule b=pi(a)) = 28, example (0, 8, 7, 6, 5, 4, 3, 2, 1)
- worst (rule b=pi(a)-1) = 22, example (0, 8, 7, 6, 5, 4, 3, 2, 1)
- worst (rule b=pi(a-1)+1) = 28, example (0, 8, 7, 6, 5, 4, 3, 2, 1)

## Search run 2026-09-16 11:09:01 (5.7 s), seed=2, iters=4000

### Reflections: fixed-a0 gap vs floor((n-1)^2/4) (exhaustive over h)
- n=10: bound=20, worst gap=0 at h=0
- n=15: bound=49, worst gap=0 at h=0
- n=20: bound=90, worst gap=0 at h=0
- n=30: bound=210, worst gap=0 at h=0
- n=40: bound=380, worst gap=0 at h=0
- n=50: bound=600, worst gap=0 at h=0
- n=60: bound=870, worst gap=0 at h=0
- n=80: bound=1560, worst gap=0 at h=0
- n=100: bound=2450, worst gap=0 at h=0

### Hill-climb (2-opt from reversed identity, 4000 moves): best gap found
- n=9: gap=1, example=(0, 7, 6, 5, 4, 3, 2, 1, 8)
- n=12: gap=3, example=(11, 2, 9, 8, 7, 6, 5, 4, 3, 1, 0, 10)
- n=15: gap=6, example=(2, 0, 13, 11, 10, 9, 8, 7, 6, 5, 4, 3, 1, 14, 12)
- n=20: gap=14, example=(4, 2, 19, 18, 16, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 3, 1, 0, 17, 15)
- n=25: gap=22, example=(6, 1, 24, 4, 22, 20, 19, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 5, 3, 2, 0, 23, 21, 18)
- n=30: gap=32, example=(9, 3, 6, 1, 27, 29, 25, 23, 22, 21, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 8, 4, 7, 5, 2, 28, 0, 26, 24, 20)
