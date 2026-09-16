/* h13_local_min.c -- does bounded-radius "joint" local optimality of the double
 * cut (q, c) certify I(pi) <= floor((n-1)^2/4) (H13-I, session 9)?
 *
 * Background (docs/notes/h13_line_model.md, session 8, Lemma C): checking
 * optimality of a cut (q, c) by shifting q alone (any k) or c alone (any k) is
 * NOT sufficient -- a counterexample at n = 7 satisfies both single-axis
 * conditions but is not the global minimum.  Session 9 asks whether a joint
 * (2D) neighborhood on the (q, c) torus fixes this: is every "king move"
 * local minimum of radius R (no cell within Linf-distance <= R has a smaller
 * inversion count) already <= floor((n-1)^2/4)?
 *
 * Reparametrization: instead of (q, c) use (q, S) with S = (q + 1 - c) mod n
 * (a bijection for fixed q; S is the value-shift applied uniformly to pi,
 * independent of q -- see line_model.py's w_j = pi(q+1+j) - (q+1-c) mod n).
 * T[q][S] = inv(w) for that cut.
 *
 * Fast recurrence (avoids the O(n^4)-per-permutation brute force of
 * line_profile.c), derived and verified against brute force by --check:
 *   inv0(S) = inversions of v_i = (pi[i] - S) mod n in natural position order
 *             (i.e. T[n-1][S], the cut q = n-1, "start at position 0").
 *   inv0(S+1) = inv0(S) + (n - 1 - 2 * posOf[S])
 *             where posOf[S] = pi^{-1}(S) (moving the value-cut by one makes
 *             the element of value S wrap from smallest to largest relabeled
 *             value; it flips relative order with all n-1 others).
 *   T[q][S]   = T[q-1][S] + (n - 1 - 2 * v[q])   (v[q] = (pi[q]-S) mod n;
 *             rotating the position-cut by one moves the front element of the
 *             current line to the back, flipping its order with all n-1
 *             others -- the same identity, applied to positions instead of
 *             values), with T[n-1][S] = inv0(S) as base case.
 * Total cost: O(n) per S, O(n) values of S -> O(n^2) per permutation (against
 * O(n^4) for the direct double loop), enabling exhaustive n = 10, 11.
 *
 * Usage: h13_local_min n [check] [R]
 *   n      -- 4 <= n <= 11 (exhaustive over all n! permutations)
 *   check  -- if the literal string "check", also verify T against brute
 *             force inversions on every cell for every permutation
 *   R      -- king-move radius for local optimality (default 2)
 * Output: one summary line to stdout:
 *   n=<n> bound=<floor((n-1)^2/4)> R=<R> worst_localmin=<v> (rank <r>) count=<n!>
 * "worst_localmin" is the maximum, over all permutations and all cells that
 * are local minima of the given radius, of the cell's inversion count.
 * Version h13_local_min-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 12

static int brute_inv(const int *pi, int n, int q, int S) {
    int w[MAXN];
    for (int j = 0; j < n; j++) {
        int v = (pi[(q + 1 + j) % n] - S) % n;
        if (v < 0) v += n;
        w[j] = v;
    }
    int inv = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (w[i] > w[j]) inv++;
    return inv;
}

static void fast_table(const int *pi, int n, int T[MAXN][MAXN]) {
    int posOf[MAXN];
    for (int i = 0; i < n; i++) posOf[pi[i]] = i;
    int inv0[MAXN];
    {
        int inv = 0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                if (pi[i] > pi[j]) inv++;
        inv0[0] = inv;
        for (int S = 0; S < n - 1; S++) {
            int p = posOf[S];
            inv0[S + 1] = inv0[S] + (n - 1 - 2 * p);
        }
    }
    for (int S = 0; S < n; S++) {
        int v[MAXN];
        for (int i = 0; i < n; i++) { int x = (pi[i] - S) % n; if (x < 0) x += n; v[i] = x; }
        T[n - 1][S] = inv0[S];
        for (int q = 0; q < n - 1; q++) {
            int prev = (q == 0) ? T[n - 1][S] : T[q - 1][S];
            T[q][S] = prev + (n - 1 - 2 * v[q]);
        }
    }
}

static int next_perm(int *a, int n) {
    int i = n - 2;
    while (i >= 0 && a[i] >= a[i + 1]) i--;
    if (i < 0) return 0;
    int j = n - 1;
    while (a[j] <= a[i]) j--;
    int t = a[i]; a[i] = a[j]; a[j] = t;
    for (int l = i + 1, r = n - 1; l < r; l++, r--) { t = a[l]; a[l] = a[r]; a[r] = t; }
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s n [check] [R]\n", argv[0]); return 2; }
    int n = atoi(argv[1]);
    if (n < 2 || n > MAXN) return 2;
    int do_check = (argc > 2 && strcmp(argv[2], "check") == 0);
    int R = (argc > 3) ? atoi(argv[3]) : 2;

    int pi[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;
    long rank = 0;
    int bound = (n - 1) * (n - 1) / 4;
    int worst = -1;
    long worst_rank = -1;
    int T[MAXN][MAXN];
    do {
        fast_table(pi, n, T);
        if (do_check) {
            for (int q = 0; q < n; q++) for (int S = 0; S < n; S++) {
                int b = brute_inv(pi, n, q, S);
                if (b != T[q][S]) {
                    printf("MISMATCH rank=%ld q=%d S=%d fast=%d brute=%d\n", rank, q, S, T[q][S], b);
                    return 1;
                }
            }
        }
        for (int q = 0; q < n; q++) for (int S = 0; S < n; S++) {
            int v = T[q][S];
            int ok = 1;
            for (int dq = -R; dq <= R && ok; dq++)
                for (int dS = -R; dS <= R; dS++) {
                    if (dq == 0 && dS == 0) continue;
                    int qq = ((q + dq) % n + n) % n, ss = ((S + dS) % n + n) % n;
                    if (T[qq][ss] < v) { ok = 0; break; }
                }
            if (ok && v > worst) { worst = v; worst_rank = rank; }
        }
        rank++;
    } while (next_perm(pi, n));

    printf("n=%d bound=%d R=%d worst_localmin=%d (rank %ld) count=%ld\n",
           n, bound, R, worst, worst_rank, rank);
    return 0;
}
