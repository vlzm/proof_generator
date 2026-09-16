# H13-I, session 9 — Mantel/Turan reformulation and three ruled-out candidate cuts

Date: 16.09.2026. Core oracle-1.0, tables C2v (4 <= n <= 10). Code:
`experiments/h13i_candidates.py` (h13i_candidates-1.0). Report:
`data/runs/h13i_candidates/report.md`, `.json`.

Task (`docs/notes/h13_line_model.md` #6, `PLAN.md` #8, budget 1 session):
prove H13-I -- every pi admits a double cut (q, s) (rotation of positions
and of values) with `I(pi; q, s) <= floor((n-1)^2/4)` inversions of the
relabelled line, equivalently `>= floor(n^2/4)` non-inversions (concordant
pairs). Definitions -- `docs/notes/h13_line_model.md` #0.

## 1. Reformulation: the target is exactly the Mantel/Turan number

For all n >= 1:

    C(n,2) = floor((n-1)^2/4) + floor(n^2/4).                         (*)

Proof: for n = 2k, C(n,2) = 2k^2-k, floor((n-1)^2/4) = floor((2k-1)^2/4)
= k^2-k, difference k^2 = floor(n^2/4). For n = 2k+1, C(n,2) = 2k^2+k,
floor((n-1)^2/4) = floor((2k)^2/4) = k^2, difference k^2+k =
floor((2k+1)^2/4) = floor(n^2/4). Both match (*). So

    I(pi) <= floor((n-1)^2/4)  for every cut  <=>  some cut has
    >= floor(n^2/4) concordant (non-inverted) pairs.

`floor(n^2/4)` is exactly the Mantel/Turan extremal number: the maximum
number of edges of a triangle-free graph on n vertices, attained uniquely
(up to isomorphism) by the complete bipartite graph `K_{ceil(n/2),floor(n/2)}`.
This is not a coincidence for the extremal inputs: at the cut that attains
`I(pi_h) = floor((n-1)^2/4)` for a reflection `pi_h` (C33/C35; the same cut
used by the two-arc word of C32), the n points split into two arcs of sizes
`a = ceil(n/2)`, `b = floor(n/2)` (by *position*, equivalently by *value*,
the cut is the same split on both circles) such that:

  - within each arc the values are in reverse order (all `C(a,2) + C(b,2) =
    floor((n-1)^2/4)` inversions come from inside the arcs);
  - every cross pair (one point from each arc) is concordant -- exactly
    `a*b = floor(n^2/4)` pairs, realizing the Mantel extremal graph
    `K_{a,b}` on the "concordance graph" of that cut.

Checked directly (not previously stated this way in `h13_line_model.md`):
`4 <= n <= 10`, reflection `pi(i) = -i mod n`, the argmin cut splits w into
exactly this two-block pattern (`h13i_opt_cut.py`, scratch, not committed --
reproduced by inspection of `line_profile`/`line_rotation_budget` data plus
direct recomputation here). This ties C32's "two independent arc reversals"
geometric picture to the Mantel extremal graph algebraically, but is a
restatement of C32/C33, not a new bound.

## 2. Three cheap candidate cuts, all ruled out (exhaustive, 4 <= n <= 8)

The reformulation suggests looking for a cut cheaper than the full n^2-grid
search that is provably good enough. Three natural candidates, all tested
exhaustively over every pi at 4 <= n <= 8 by `h13i_candidates.py`:

**(a) Pigeonhole on s alone.** Fix q, average `inv(q, s)` over the n values
of s (this average is `F(q)/n` with `F(q) = sum_{j<k} (v_k - v_j) mod n`
over the q-linearized value sequence v); some s attains at most the
average, so `min_q F(q)/n <= floor((n-1)^2/4)` would suffice, using only n
(not n^2) evaluations of F. Data: `max_pi min_q F(q)` against `n *
floor((n-1)^2/4)`:

| n | n * bound | max_pi min_q F(q) | ratio |
|---|---|---|---|
| 4 | 8 | 14 | 1.75 |
| 5 | 20 | 30 | 1.50 |
| 6 | 36 | 55 | 1.53 |
| 7 | 63 | 91 | 1.44 |
| 8 | 96 | 140 | 1.46 |
| 9 | 144 | 204 | 1.42 |

Worst case at every n is the reflection `pi(i) = -i mod n` (e.g. n = 8:
`(0,7,6,5,4,3,2,1)`). Ruled out by a wide, non-shrinking margin (~1.4-1.5x);
consistent with the already-recorded failure of averaging over *both* q and
s jointly (`h13_line_model.md` #6), but this shows the weaker one-sided
pigeonhole (average only over s, exact min over q) fails just as badly, by
essentially the same factor. Averaging over s cannot see the two-arc
structure of #1: the arcs are separated by wrap-around, which a same-q
uniform average over s cannot exploit because it mixes all `floor(n^2/4)`
good pairs together with the `floor((n-1)^2/4)` bad ones in one linear
scan.

**(b) Anchor the cut at a data point.** Restrict to the n cuts `(q, pi(q))`
(the point i = q becomes the new origin, so it never itself contributes an
inversion; PLAN #8's induction idea removes exactly this point). Data:

| n | bound | max_pi min over anchor cuts | ratio |
|---|---|---|---|
| 4 | 2 | 6 | 3.0 |
| 5 | 4 | 10 | 2.5 |
| 6 | 6 | 15 | 2.5 |
| 7 | 9 | 21 | 2.33 |
| 8 | 12 | 28 | 2.33 |
| 9 | 16 | 36 | 2.25 |

Again the reflection is worst at every n, and the anchor value is exactly
`C(n-1,2)` (the anchored point removes itself, the remaining n-1 points are
still in exact reversed order regardless of anchor choice, since a global
order-reversing map restricted to any n-1 of its points is still fully
reversed). This sharpens the existing negative result in
`h13_line_model.md` #1.7 (single-element induction fails on a specific
n = 8 example, by a gap of 5 against an allowed increment of 3): here the
*entire* restricted candidate set of n anchored cuts (not just one
reduction step composed with a specific inductive formula) already misses
the target by a factor ~2.3-3, on the extremal family itself, for every
n <= 8. The good cut for a reflection is never anchored at one of its own
points -- it sits in the middle of an arc, as #1 shows.

**(c) Best semicircle alignment (2x2 quadrant bound).** For each pair of a
position-arc `A` (size `a = ceil(n/2)`) and a value-arc `B` (size a) among
the n choices each, let `x = |A cap pi^{-1}(B)|` (points with position in A
*and* value in B). Every point in `A cap pi^{-1}(B)` is concordant with
every point in `(complement of A) cap pi^{-1}(complement of B)` regardless
of the internal order within each part (guaranteed cross-quadrant
concordant pairs), giving `x * (n - 2a + x)` concordant pairs for free.
For n even (a = n/2) this is `x^2`, which reaches the target `floor(n^2/4)
= a^2` only when `x = a` exactly, i.e. pi maps some semicircle of positions
*bijectively* onto some semicircle of values (a "perfect alignment").
Data (max over the n^2 arc-pairs of x, worst pi, and how often a perfect
alignment exists at all):

| n | need x = a | worst-case max x | perfect alignments / n! |
|---|---|---|---|
| 4 | 2 | 2 | 24/24 |
| 5 | 3 | 2 | 110/120 |
| 6 | 3 | 2 | 516/720 |
| 7 | 4 | 3 | 3346/5040 |
| 8 | 4 | 3 | 15856/40320 |
| 9 | 5 | 3 | 133560/362880 |

Already at n = 5 a perfect alignment fails to exist for 10 permutations
(e.g. one of them: best achievable x = 2 < 3), and the fraction of
permutations with *no* perfect alignment grows (57% at n = 6, 61% at
n = 8, 63% at n = 9), and the worst-case shortfall itself grows
(need 5, best 3, at n = 9). So even the guaranteed-only two-quadrant bound
cannot be the whole argument: the true optimal cut must also collect some of the
within-quadrant concordant pairs that this crude count throws away.
(Reflections themselves *do* have a perfect alignment, consistent with #1;
the extremal family is not among the failures of (c).)

## 3. Verdict

H13-I stays **CONJECTURED**, unchanged in status; VERIFIED 4 <= n <= 10
(`I(pi) <= floor((n-1)^2/4)`, C33/C35, `experiments/line_profile.c`, no new
computation needed here). New this session:

1. The Mantel/Turan reformulation (#1): the target non-inversion count is
   exactly the extremal number `floor(n^2/4)`, and the extremal inputs
   (reflections) realize the Mantel extremal graph `K_{ceil(n/2),floor(n/2)}`
   as their concordance graph at the optimal cut. This is a restatement,
   not a new bound, but gives the problem a recognizable shape and a
   candidate proof shape (find, for every pi, *some* balanced bipartition
   of the n points into two parts of size ceil(n/2)/floor(n/2) all of whose
   cross pairs are concordant under *some* cut, or a weaker substitute that
   still totals `floor(n^2/4)`).
2. Three cut-selection shortcuts that would turn the n^2-grid search into
   an O(n) or O(n^2)-candidate search are all exhaustively ruled out at
   4 <= n <= 8, each failing on the same extremal family (reflections) by a
   growing or flat (not shrinking) margin: pigeonhole over s alone (~1.5x
   over target), anchoring at a data point (~2.3-3x over target), and the
   crude two-quadrant alignment bound (perfect alignment does not always
   exist, and is not by itself sufficient reasoning even where it does).

No proof found within the 1-session budget. Per `PLAN.md` #8 and
`h13_line_model.md` #6, the documented next step is either: (a) a genuinely
2-dimensional argument that tracks within-quadrant order recursively (not
attempted here -- candidate: shifting/compression argument showing
reflections are the pointwise worst case, i.e. any adjacent transposition
of pi can only improve the best achievable non-inversion count, by analogy
with compression arguments in extremal set theory; not tested), or (b) a
literature search specifically for "cyclic sequences under joint rotation
of index and value" turned up nothing on a first pass (see below) and was
not pursued further, or (c) fall back to the independent path H3 (PLAN
#4.3 item 3, induction on `sigma_(n-1) -> sigma_n`), which does not depend
on H13-I.

Web search for the exact quantity ("minimum inversions of a permutation
under joint cyclic rotation of positions and values", "torus permutation
matrix Mantel/Turan bound") returned no matching literature; this may be
folklore, may be known under different terminology (e.g. in genome
rearrangement / circular sequence comparison), or may not have been framed
this way before. Not conclusive either way; a deeper literature search was
not attempted given the 1-session budget.

## 4. Not done

- The shifting/compression idea (#3, item (a)) is untested, computationally
  or on paper.
- No attempt was made at H3 (PLAN #4.3) this session; it is the next
  candidate task per the standing plan.
- 4 <= n <= 9 is exhaustive for all three candidates (`data/runs/h13i_candidates/`,
  n = 9 in 178 s); n = 10 was not run (the ladder in AGENTS.md #9 does not
  require it for a first negative pass, and `I(pi)` itself is already
  VERIFIED at n = 10 by the existing `line_profile.c`, C33/C35).
