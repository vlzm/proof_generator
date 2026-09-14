# H3 induction reduction — report (sessions 10-11, 14.09.2026)

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
| 8 | 7 | 9 | 8 | **EXCEEDS (+1)** |
| 9 | 8 | 10 | 8 | ok |
| 10 | 9 | 12 | 10 | **EXCEEDS (+1)** |
| 11 | 10 | 12 | 11 | **EXCEEDS (+1)** |
| 12 | 11 | 15 | 13 | **EXCEEDS (+2)** |

Logs: `reduction_n5_10.log`, `reduction_n11.log`, `reduction_n12.log`. The
n=12 run (44 CPU-min, O(n^4) per permutation) was launched in the
background and not originally awaited when this report's verdict was
first written; it finished afterward and is folded in here as a same-day
addendum. It matters: the excess is not flat at +1, it grows to +2 at
n=12.

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
11, 12 (not an isolated small-n artifact), and the excess is not flat —
it doubles from +1 (n=8,10,11) to +2 (n=12). This mirrors the earlier H10/H11
finding (independent per-cycle processing loses Theta(n)) and the H13-I
finding (single-coordinate/averaged cut selection is insufficient): every
attempt so far to decompose the problem into independent smaller pieces
loses exactly where O(1) precision is needed. Session 11 tested the
two-level induction step (remove a PAIR of elements at once, see below);
it too was refuted, closing H3's 2-session budget.

## Session 11: two-element removal (k=2)

`experiments/h3_reduction_pair.c` (h3_reduction_pair-1.0): remove a pair
of values {v1,v2} at once (C(n,2) choices, no rotation), compare
d_n(pi) to d_{n-2}(reduced); target is the two-step telescoped budget
2n-3 = (n-1)+(n-2).

| n | target = 2n-3 | worst | verdict |
|---|---|---|---|
| 6 | 9 | 9 | ok (exact) |
| 7 | 11 | 11 | ok (exact) |
| 8 | 13 | 13 | ok (exact) |
| 9 | 15 | 15 | ok (exact) |
| 10 | 17 | 17 | ok (exact) |
| 11 | 19 | 19 | ok (exact) |
| 12 | 21 | 22 | **EXCEEDS (+1)** |

Logs: `pair_n6_10.log`, `pair_n11.log`, `pair_n12.log`. Notably, 6 <= n <=
11 hit the target EXACTLY (not just below) — a much cleaner run than k=1
ever produced — before breaking at n=12 with counterexample
`pi=(4,3,6,5,0,11,2,1,8,7,10,9)` (cycle type: two 2-cycles, one 6-cycle,
two fixed points; worst=22>21).

## Final conclusion (H3 closed)

Both k=1 (session 10) and k=2 (session 11) "delete k elements, best
choice, no rotation cost" reductions fail eventually: k=1 at n=8, k=2 at
n=12. Larger k buys more room before the first counterexample, but this
is not a path to C15 with a *fixed* k — the proof needs to work at every
n, and letting k grow with n degenerates the induction into the original
unsolved problem. Registry: C38 (k=1), C39 (k=2). H3's 2-session budget
is spent (AGENTS.md rule 12); next front is auditing the S1 lower bound
(C8, PLAN S3.2), per `docs/notes/h3_induction.md` S7.
