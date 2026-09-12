# A carrier-count/gap improvement of the LRX diameter bound

Date: 2026-09-11. Author: GitHub Copilot, with an independent mathematical
subagent audit. This is a new derivation relative to the supplied project
proof, not a claim of priority over all literature. No Telegram collaboration
was used, per the user's later instruction.

## Result and dependencies

Use the generators, left-to-right word convention, and notation of the
[supplied proof](lrx_2n3_proof_checked.md). Put

$$
P=\lfloor n^2/4\rfloor,\qquad
A_n=2P+n-\frac{2P+11n-11}{3n}.
$$

**Theorem.** For every integer $n\ge4$,

$$
\boxed{D_n\le\lfloor A_n\rfloor-1
=2P+n-1-\left\lceil\frac{2P+11n-11}{3n}\right\rceil.}
\tag{T1}
$$

In particular,

$$
D_n\le\binom n2+\left\lfloor\frac{4n}{3}\right\rfloor-4.
$$

A second bound, useful at larger orders, is

$$
\boxed{D_n\le\left\lfloor A_n-
\max\left\{1,\,3-\frac{138}{n}\right\}\right\rfloor.}
\tag{T2}
$$

Thus the supplied integer bound improves by at least two for $n\ge138$,
and by at least three for $n\ge852$ (rounding justification below).
These are constant improvements, **not a proof of the exact diameter**.

The argument uses the proved per-cycle word estimate, commuting-cycle
assembly, and averaging identities from Sections 1–5 of the supplied proof.
Those have been read and checked algebraically, and the cycle words have
been independently executed in the finite tests below. No lower-bound
literature is needed for this upper bound.

## 1. Slack on short cycles

For a nontrivial cycle $C$ of length $k$, let $F,M,t,E=k-2t$ have the
meanings in the supplied proof. Its internal word has length at most $2F-S$,
where

$$
S=2M+E-1,\qquad F+k\le nM+E,\qquad
B:=S-\frac{5(F+k)}{3n}\ge0.
\tag{1}
$$

**Low-load lemma.** If $M\le2$, then $S\ge k+1$. Consequently, always

$$S\ge\min\{k+1,\,5+(k\bmod2)\}.\tag{2}$$

**Proof.** Adjacent edge loads differ by $0$ or $\pm2$, so all loads
have the same parity. If $M=1$, all loads are one and there are no sign
changes; $E=k$, giving $S=k+1$.

If $M=2$ and all signs agree, $E=k$ and the conclusion is stronger.
Otherwise the loads are zero or two and a sign change forces a zero.
The support of the paths is connected because their concatenation is a
single closed walk. Hence it is one proper arc of the circle. There is
exactly one transition from load zero to two. The load-difference formula
in the supplied proof identifies every $-\to+$ change with such a
transition. Thus $t=1$, $E=k-2$, and $S=k+1$.
For $M\ge3$, $S\ge5+E\ge5+(k\bmod2)$. ∎

**Short-cycle lemma.** For $n\ge12$ and $2\le k\le7$,

$$
B\ge b_n:=\min\left\{\frac53-\frac{10}{n},\,
1-\frac5{3n}\right\}.
\tag{3}
$$

Here $b_{12}=5/6$, and $b_n=1-5/(3n)$ for $n\ge13$.

**Proof.** If $k\le5$, (2) and $F\le kn/2$ imply

$$B\ge1+\frac{k(n-10)}{6n}\ge1.$$

For $k=6$, unless $M=3,E=0$, (2) and parity imply $S\ge7$,
and therefore $B\ge2-10/n$. In the exceptional case there are three
positive and three negative arrows. The clockwise load minus the
counterclockwise load on every edge equals one constant integer $q$:
this follows from flow conservation for a closed walk. Moreover
$q\equiv M\pmod2$ and $|q|\le M=3$. The values $q=\pm3$ would force
all edge traversals into one direction, contrary to the mixed signs.
So $q=\pm1$. If $U,V$ are the total positive and negative lengths,
$|U-V|=n$ and $U,V\le3n/2$. Hence $F=U+V\le2n$, yielding
$B\ge5/3-10/n$.

For $k=7$, either $S\ge8$, giving
$B\ge13/6-35/(3n)\ge1$ for $n\ge12$, or $S=6$.
The latter forces $M=3,E=1$ by (2), so $F+k\le3n+1$ and
$B\ge1-5/(3n)$. ∎

For a permutation with $K$ nontrivial cycles, let $J$ count cycles of
length at most seven, and $N$ be the number of nonfixed vertices. Then

$$
n\ge N\ge2J+8(K-J),\qquad
\sum_C B(C)\ge b_n\max\left\{0,
\left\lceil\frac{8K-n}{6}\right\rceil\right\}.
\tag{4}
$$

## 2. Positioning walks can omit gaps

Fix a shift $c\ne0$, let $d=\delta_n(0,c)$, and choose one of the
permitted carriers for each cycle of $f_c$. Their positions are distinct.
In a chosen shortest $0$-to-$c$ arc, insert these carriers and the two
endpoints as required vertices. The largest gap $g$ between consecutive
required vertices has length

$$g\ge\left\lceil\frac d{K+1}\right\rceil.$$

A walk can visit all carriers while omitting that gap. If its endpoints
have coordinates $x<y$ along the chosen short arc, travel $0\to x\to0$,
then along the complementary arc to $c$, then $c\to y\to c$. Its length is

$$2x+(n-d)+2(d-y)=n+d-2g.$$

Perform each cycle at its carrier's first visit; the disjoint cycles
commute, so this is a valid assembly of the supplied internal words.
Compared with $H(c)=n+d-2$, the guaranteed external saving is

$$
2(g-1)\ge2\left\lfloor\frac{d-1}{K+1}\right\rfloor.
\tag{5}
$$

The same construction on the other arc gives an alternative route of
length $2n-d-2g_\ell$, where $g_\ell$ is its largest gap. At $c=0$
we retain the original full-circuit bound $H(0)=n$.

Two elementary improvements will also be useful at small orders:

- For $K=0$, the actual route has length $d$, saving $H(c)-d$.
- For $K=1$ and $c\ne0$, the route length is at most $n-d$, saving
  at least $2d-2$. If the carrier lies on the short arc, use that arc.
  If it lies on the long arc, the sum of the two shortest distances
  from the endpoints to the carrier is at most the long arc's length.
- At an antipodal endpoint, if $K<n-2$, some interior vertex of the two
  arcs is not a carrier. One arc has a gap of length at least two, so
  the better of the two gap routes saves at least two.

## 3. Averaging the combined savings

Let $R_c$ be total saving from the exact cycle slacks $\sum B$ and from
an improved positioning route, relative to the original estimate for
shift $c$. The nonnegativity in (1) and the original route give $R_c\ge0$.
Every permutation can be sorted with cost at most

$$
2F_c-\frac5{3n}(F_c+N_c)+H(c)-R_c.
$$

The supplied identities $\operatorname{avg}F_c=P$,
$\operatorname{avg}N_c=n-1$, and
$\operatorname{avg}H(c)=n-2+(P+2)/n$ imply

$$D_n\le\left\lfloor A_n-\frac1n\sum_cR_c\right\rfloor\tag{6}$$

whenever the bound on the average saving holds uniformly over inputs.
There is no assumption that the cycle counts of different shifts are
independent or identical.

For $n\ge12$, put $q_0=\lfloor(n+6)/8\rfloor$. At every endpoint
with $d\ge q_0+2$:

- if $K\le q_0$, (5) saves at least two;
- if $K>q_0$, (4) gives at least two short cycles, saving at least $2b_n$.

As $b_n\le1$, both cases guarantee $R_c\ge2b_n$. There are exactly
$n-2q_0-3$ such endpoints. Thus

$$
\frac1n\sum_cR_c\ge G_n:=
\frac{2b_n}{n}\left(n-2\left\lfloor\frac{n+6}{8}\right\rfloor-3\right).
\tag{7}
$$

For $n\ge22$,

$$
G_n\ge\left(1-\frac5{3n}\right)\left(\frac32-\frac9n\right)
=\frac32-\frac{23}{2n}+\frac{15}{n^2}\ge1.
$$

The last expression is increasing for $n\ge22$ and equals $122/121$
at $n=22$. Directly, $G_{20}=121/120$ and $G_{21}=464/441$.
This proves (T1) for all $n\ge20$.

## 4. Small orders without relying on a diameter table

For even $n=2m$ and $k\ge3$, antipodal arrows in a cycle form a
matching: two of them sharing a vertex would force a two-cycle.
There are at most $\lfloor k/2\rfloor$ such arrows; all others have
length at most $m-1$. Consequently

$$F\le k(m-1)+\lfloor k/2\rfloor.\tag{8}$$

Here are uniform lower bounds on the average in (6):

| $n$ | Average saving at least | $n$ | Average saving at least |
|---:|---:|---:|---:|
| 4 | $5/8$ | 12 | $25/36$ |
| 5 | $1$ | 13 | $136/169$ |
| 6 | $10/9$ | 14 | $37/42$ |
| 7 | $16/21$ | 15 | $128/135$ |
| 8 | $13/12$ | 16 | $129/128$ |
| 9 | $62/81$ | 17 | $920/867$ |
| 10 | $7/10$ | 18 | $49/54$ |
| 11 | $104/121$ | 19 | $1040/1083$ |

Entries for $12\le n\le19$ are (7). For the others:

- **$n=4$.** Equations (2), (8) give $B\ge1/2$ per cycle. At the
  antipodal endpoint, one carrier saves two externally, while two cycles
  give slack at least one. Other endpoints have saving at least $1/2$.
  Average: at least $(3/2+1)/4=5/8$.
- **$n=5$.** $F+k\le3k$, $S\ge k+1$, so $B\ge1$ per cycle.
- **$n=6$.** With $K\ge2$, each cycle has $k\le4$ and
  $B\ge7/9$ by (2), (8). Thus at $d=2$ the saving is at least $14/9$.
  At $d=3$ add two for the antipodal gap, giving $32/9$.
  If $K=1$, external savings $2$ and $4$ already suffice.
  Average at least $(2(14/9)+32/9)/6=10/9$.
- **$n=7$.** For $k\le5$, (2) gives $B\ge23/21$.
  At $k=6,M=3,E=0$, three positive steps of length at most three
  and three negative steps of length at least one cannot have
  $|U-V|=7$, so this exceptional case is impossible. Otherwise
  $k=6$ gives $B\ge9/7$. At $k=7$, the only delicate case is
  $S=6,M=3,E=1$, giving $B\ge16/21$ by (1).
- **$n=8$.** With $K\ge2$, total slack is at least $4/3$.
  For $k=2,3,4,5$, respective lower bounds are
  $11/12,31/24,5/4,17/12$, from (2), (8).
  For $k=6$ the exceptional case gives $5/12$ using $F\le16$;
  all other cases give at least $11/8$. The smallest possible
  combined bound is therefore from $6+2$, namely $4/3$.
  At $d=2,3$ use this; at $d=4$ add two for the antipodal gap.
  If $K=1$, the external savings suffice instead.
  Average at least $(4(4/3)+10/3)/8=13/12$.
- **$n=9,10,11$.** If $K\ge2$, some cycle has $k\le5$.
  From (2) and $F\le k\lfloor n/2\rfloor$, its slack is at
  least $31/27,1,13/11$, respectively. If $K=1$ and $d\ge2$,
  the external saving is at least two instead. There are $n-3$
  endpoints with $d\ge2$, giving the entries shown.

In every case $K=0$, the exact route of length $d$ gives at least the
saving needed in the relevant argument. Each table entry less than one
is **strictly greater than the fractional part of $A_n$**. Equation (6)
therefore reduces $\lfloor A_n\rfloor$ by at least one, proving (T1).

## 5. Asymptotically three saved moves

For $n\ge13$, put $Q_n=\lfloor(n+36)/8\rfloor$ and define

$$
T_n(d)=2\left\lfloor\frac{d-1}{Q_n+1}\right\rfloor\quad(d\ge1),
\qquad T_n(0)=0.
$$

If $K>Q_n$, (4) forces $J\ge7$, so the cycle slack alone is at
least $7b_n\ge6$. If $K\le Q_n$, (5) guarantees $R_c\ge T_n(d)$.
The surrogate $T_n(d)$, for $d\le\lfloor n/2\rfloor$, is at most six;
hence in both cases $R_c\ge T_n(d)$. (Actual savings can be larger.)
Counting the three thresholds gives

$$
\frac1n\sum_cR_c\ge Z_n:=\frac2n
\sum_{j=1}^3\bigl(n-2j(Q_n+1)-1\bigr)_+.
\tag{9}
$$

Using $(x)_+\ge x$ and $Q_n\le(n+36)/8$,

$$Z_n\ge6-\frac{24Q_n+30}{n}\ge3-\frac{138}{n}.$$

Together with (T1), this proves (T2); for $4\le n\le12$ its maximum
is just one. The exact bounds (7), (9) can also be used directly.

For $n\ge852$, the fractional part of $A_n$ is at most
$5/6+23/(6n)$. To check this, use the parity expansions in Section 5
of the supplied proof and the three residues of $n$ modulo three;
for even $n$ the largest fraction is only $2/3+11/(3n)$.
Thus

$$\{A_n\}+\frac{138}{n}\le\frac56+\frac{851}{6n}<1.$$

It follows that $\lfloor A_n-3+138/n\rfloor\le\lfloor A_n\rfloor-3$.

## 6. Verification and limitations

- A separate mathematical subagent audited Sections 1–3 and 5. It found
  no substantive error and requested the clarifications that the gains
  are lower bounds and the six-move cap concerns the surrogate only.
- [verify_lrx_gap_bound.py](verify_lrx_gap_bound.py) independently
  constructs and executes internal cycle words, including all shortest-path
  tie choices and all permitted carriers for $4\le n\le8$.
- It compares the gap-route formula with BFS on the circle plus a visited
  carrier mask, for every carrier subset of size at most $n/2$ and every
  endpoint, for $4\le n\le10$.
- Separate exact-rational finite certificates check the small orders using
  relaxed cycle parameters and all partitions. These are stronger than
  parts of the simple table above, but are not needed for its proof.
- Full Cayley BFS independently confirms $D_n=\binom n2$ for $4\le n\le8$.
  This is finite evidence only. Experiment counts and commands are in
  [LOGS.md](LOGS.md).

**What remains open here.** Neither equality $D_n=\binom n2$ for every
$n\ge4$ nor an improvement of the linear-in-$n$ excess has been proved.
The present inequalities allow zero charged slack on cycles of length
eight or more and $K$ of order $n$, explaining why the guaranteed gap
saving is only constant. A larger improvement likely needs a new cycle
construction or constraints coupling different shifts $f_c$.