"""H13-I candidate: shift-recursion structure on the double cut (q, c) and the
"concordant pairs" reformulation.

Recall (docs/notes/h13_line_model.md): a double cut (q, c) turns pi into the
line w_j = pi(q+1+j) - (q+1-c) mod n, j = 0..n-1.  H13-I claims
min_{q,c} inv(w) <= floor((n-1)^2/4) for every pi.

This module checks two things found while looking for a proof of H13-I
(session 9, no proof found; see docs/notes/h13_line_model.md §7):

1. Claim A (fixed position cut suffices): is min_c inv(w) <= floor((n-1)^2/4)
   already true for q fixed (e.g. q = n-1, i.e. w built directly from pi
   without rotating positions)?  REFUTED: worst case at n = 7 is 10 > 9,
   attained at pi = (0,5,4,3,2,1,6) -- the same permutation that refutes the
   local-optimality conditions of Lemma C (h13_line_model.md §4).  So fixing
   q is no easier than the general problem; genuine use of both q and c is
   needed, confirmed for 4 <= n <= 9 (n = 9 exhaustive is slow in Python,
   sampled here; C33/C35 already establish n <= 10 for the full min_{q,c}).

2. Clean shift recursion.  Reparametrising the value cut by an INDEPENDENT
   anchor c' (c' = q+1-c, so {(q,c)} and {(q,c')} range over the same n^2
   cuts for fixed q) decouples the two moves:
     Op_Q (q -> q+1, c' fixed):  inv changes by (n-1) - 2*w_0
        (w_0 = current value at line position 0; this is a pure "move the
        front element to the back" step, no relabelling).
     Op_C (c' -> c'+1, q fixed): inv changes by (n-1) - 2*pos_0
        (pos_0 = position of the element with value 0; a pure "relabel all
        values down by 1 mod n" step).
   Both identities are exact (proved by direct inversion counting: an
   element with value m at the front contributes m inversions before the
   move and n-1-m after; symmetric for value 0 under relabelling) and
   verified below on random instances.  This gives inv(q,c') as a height
   function on the n x n torus with a prescribed multiset of row/column
   steps {(n-1)-2k : k=0..n-1}, but a coordinate-descent argument on this
   function cannot alone give H13-I: local minima can exceed the target
   (Lemma C's n=7 counterexample is already a coordinate-wise local min).

3. Concordant-pairs reformulation.  Since inv(w) + noninv(w) = C(n,2), H13-I
   is equivalent to: max_{q,c} noninv(w) >= C(n,2) - floor((n-1)^2/4), and
   the right-hand side equals floor(n^2/4) exactly for every n >= 1 (checked
   below).  So H13-I <=> "some cut makes at least floor(n^2/4) = B_n pairs
   concordant".  Per-pair count: for an unordered pair {i,j} with
   a = (j-i) mod n, b = (pi(j)-pi(i)) mod n (a, b in 1..n-1, using the
   representative with a in 1..n-1), the number of the n^2 cuts making the
   pair concordant is g(a,b) = (n-a)(n-b) + a*b.  Averaging over all n^2
   cuts fails exactly where H13-I is tight: on reflections (where
   I(pi_h) = floor((n-1)^2/4) exactly, C33/C35), the average of inv(w) over
   all cuts is close to n^2/3 -- well above the target -- so
   min_{q,c} inv <= average does not apply to the very permutations that
   need it; a second-moment argument, a non-uniform distribution over cuts,
   or an explicit construction is still needed.

No claim in this file proves H13-I; it is a record of the two structural
facts (recursion identities, concordant-pairs dual form) for the next
session, plus the exhaustive refutation of Claim A.

Usage: python3 experiments/toric_shift_recursion.py
"""

import itertools
import random

VERSION = "toric_shift_recursion-1.0"


def inv(seq):
    n = len(seq)
    c = 0
    for i in range(n):
        for j in range(i + 1, n):
            if seq[i] > seq[j]:
                c += 1
    return c


def target(n):
    return ((n - 1) ** 2) // 4


def line(pi, n, q, cp):
    return [(pi[(q + 1 + j) % n] - cp) % n for j in range(n)]


def verify_recursions(trials=5000, seed=1, nmax=14):
    rng = random.Random(seed)
    for _ in range(trials):
        n = rng.randint(3, nmax)
        pi = list(range(n))
        rng.shuffle(pi)
        q = rng.randint(0, n - 1)
        cp = rng.randint(0, n - 1)
        seq = line(pi, n, q, cp)
        i0 = inv(seq)

        seq_q1 = line(pi, n, (q + 1) % n, cp)
        pred_q1 = i0 + (n - 1) - 2 * seq[0]
        assert inv(seq_q1) == pred_q1, ("Op_Q", n, q, cp)

        pos0 = seq.index(0)
        seq_c1 = line(pi, n, q, (cp + 1) % n)
        pred_c1 = i0 + (n - 1) - 2 * pos0
        assert inv(seq_c1) == pred_c1, ("Op_C", n, q, cp)
    print(f"[recursions] {trials} random (n,q,c') trials, n<=({nmax}): both identities hold exactly.")


def refute_claim_a(nmax=8):
    print("[claim A] fixed q (q = n-1, i.e. w built directly from pi), min over c only:")
    worst_examples = {}
    for n in range(3, nmax + 1):
        worst = -1
        worst_pi = None
        for pi in itertools.permutations(range(n)):
            best_c = min(inv(line(list(pi), n, n - 1, c)) for c in range(n))
            if best_c > worst:
                worst, worst_pi = best_c, pi
        worst_examples[n] = (worst, worst_pi)
        status = "OK" if worst <= target(n) else "EXCEEDS target -> Claim A false"
        print(f"  n={n}: target={target(n)}, worst min_c inv={worst} at {worst_pi}  [{status}]")
    return worst_examples


def verify_dual_target(nmax=200):
    ok = True
    for n in range(1, nmax + 1):
        total_pairs = n * (n - 1) // 2
        dual = total_pairs - target(n)
        if dual != (n * n) // 4:
            ok = False
            print(f"  MISMATCH at n={n}: C(n,2)-target={dual}, floor(n^2/4)={(n*n)//4}")
    print(f"[dual form] C(n,2) - floor((n-1)^2/4) == floor(n^2/4) for all 1 <= n <= {nmax}: {ok}")


def average_fails_on_reflections(nmax=20):
    """Averaging over all n^2 cuts fails exactly where H13-I is tight: on
    reflections pi_h(i) = h - i mod n (equality I(pi_h) = floor((n-1)^2/4),
    C33/C35), the average of inv(w) over all n^2 cuts is close to n^2/3,
    well above the target -- so min_{q,c} inv <= average does not apply to
    the very permutations that need it.  Checked directly (inv form) and
    in the dual concordant-pairs form (same conclusion, since the two are
    complementary: avg_inv + avg_concordant = C(n,2))."""
    print("[dual averaging] reflections pi_h (h=0): avg_inv over all n^2 cuts vs target:")
    for n in range(4, nmax + 1):
        pi = [(-i) % n for i in range(n)]
        total = 0
        mn = 10 ** 9
        for q in range(n):
            for cp in range(n):
                v = inv(line(pi, n, q, cp))
                total += v
                mn = min(mn, v)
        avg = total / (n * n)
        t = target(n)
        print(f"  n={n}: avg_inv={avg:.2f}, target={t}, true min_inv={mn} "
              f"[{'averaging would fail here' if avg > t else 'averaging OK here'}]")


if __name__ == "__main__":
    print(f"{VERSION}\n")
    verify_recursions()
    print()
    refute_claim_a(nmax=8)
    print()
    verify_dual_target()
    print()
    average_fails_on_reflections()
