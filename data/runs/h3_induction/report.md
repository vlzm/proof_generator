# H3 induction reduction — report (session 10, 14.09.2026)

Goal: test the H3 induction hypothesis (PLAN.md S4.3.3 / S8) — a
"delete one array element, relabel" reduction of pi in S_n to pi' in
S_{n-1} satisfying `d_n(pi) <= d_{n-1}(pi') + (n-1)` for the best choice
of element removed, which would telescope (B_n - B_{n-1} = n-1 always,
base case D_4 = B_4 = 6) to a zero-slack proof of D_n <= B_n (C15).
Full writeup: `docs/notes/h3_induction.md`. Verdict: REFUTED for this
reduction family — counterexamples at n = 8, 10, 11.

## Command / versions

- `gcc -O2 -o h3_reduction experiments/h3_reduction.c` (h3_reduction-1.0),
  core oracle-1.0.
- Certified tables `data/tables/dist_n{4..10}.bin` (in git); n = 11, 12
  reproduced locally via `gcc -O2 -o bfs_fast exact/bfs_fast.c &&
  ./bfs_fast 11 dist_n11.bin && ./bfs_fast 12 dist_n12.bin` (16 s / 4 min
  on this machine) — not stored in git per existing convention
  (`docs/notes/repo_state.md`), reproducible and only needed for this
  one-off exhaustive check.
- `./h3_reduction n tables_dir` for n = 5..11 (12 launched, background).

## Method

For every pi in S_n (exhaustive, all n!), for every choice of element
value v to remove (n choices — delete the array entry equal to v, close
the gap, relabel remaining values down by 1 to land in S_{n-1}), compute
`extra = d_n(pi) - d_{n-1}(reduced)`. Two variants:
- **no_rotation**: use the reduced array exactly as produced by deletion
  (this is what an actual algorithm would have to pay for, since any
  further repositioning of the head costs real L/R moves in the n-sized
  problem too).
- **with_free_rotation**: additionally minimize over all n-2 rotations of
  the reduced array (an optimistic upper bound — charges nothing for the
  repositioning), to test whether the reduction family could *ever* work
  even in the most generous accounting.

Report `worst_X = max_pi min_{v[,rotation]} extra` for each variant; the
H3 recurrence needs `worst_with_free_rotation <= n-1` (a necessary
condition; the no_rotation number is what an actual construction would
realistically have to beat).

## Results

| n | target = n-1 | worst_no_rotation | worst_with_free_rotation | verdict |
|---|---|---|---|---|
| 5 | 4 | 4 | 4 | ok |
| 6 | 5 | 5 | 5 | ok |
| 7 | 6 | 7 | 6 | ok |
| 8 | 7 | 9 | 8 | **EXCEEDS** |
| 9 | 8 | 10 | 8 | ok |
| 10 | 9 | 12 | 10 | **EXCEEDS** |
| 11 | 10 | 12 | 11 | **EXCEEDS** |
| 12 | 11 | — | — | launched in background (~38 CPU-min estimated, O(n^4) per permutation), not awaited -- not needed, three counterexamples (n=8,10,11) already settle the verdict |

Logs: `reduction_n5_10.log`, `reduction_n11.log`. Reproduce n = 12 with
`./bfs_fast 12 dist_n12.bin && ./h3_reduction 12 <dir with n=11,12 tables>`
if a future session wants that extra data point.

Smallest counterexample (n, then lexicographic rank — AGENTS.md rule 6):
n = 8, pi = (3,6,0,7,4,5,2,1), d_8(pi) = 23; best achievable via any
element removal + free rotation in the n=7 subproblem is d_7 >= 15, i.e.
23 - 15 = 8 > 7 = target. The sigma_8-conjugation automorphism
(PROBLEM S4.3) does not help: pi is a fixed point of it.

## Conclusion

The "delete one array element" reduction family (any choice of element,
any rotation of the result, rotation cost uncharged — the most generous
version) cannot satisfy the zero-slack recurrence needed for a clean
induction proof of C15: it fails starting at n = 8 and again at n = 10,
11 (not an isolated small-n artifact). This mirrors the earlier H10/H11
finding (independent per-cycle processing loses Theta(n)) and the H13-I
finding (single-coordinate/averaged cut selection is insufficient): every
attempt so far to decompose the problem into independent smaller pieces
loses exactly where O(1) precision is needed. Next steps (not attempted
this session): remove a whole cycle instead of one element (same
independence risk as H11), a two-level induction step, a non-zero-slack
recurrence (informative but not a C15 proof by itself), or pivot to the
S1/C8 lower-bound audit as an alternative front.
