/* h13i_reduction.c — exhaustive check of the "remove k elements" induction
 * step for H13-I (docs/notes/h13_line_model.md §6): does every pi in S_n have
 * a k-subset K of positions such that
 *     I(pi) - I(contract(pi, K)) <= f(n) - f(n-k),   f(m) = floor((m-1)^2/4),
 * i.e. the exact per-step increment needed to prove I(pi) <= f(n) by
 * induction on n (step k)?  I(pi) = min over the n^2 double cuts (q, c) of
 * the inversion count of the relabelled line (see experiments/line_model.py).
 *
 * Session 8 showed k=1 fails at n=8 (docs/notes/h13_line_model.md §1.7,
 * pi = (0,5,2,7,4,1,6,3)); this program checks k=1,2,3 uniformly and
 * exhaustively for 4 <= n <= NMAX, with the correct (tight) budget
 * f(n)-f(n-k) rather than an ad hoc one.
 *
 * Method: precompute I(rho) for every rho in S_{n-k} (same O(n^2) cuts x
 * O(n^2) inversions per permutation as line_profile.c), indexed by Lehmer
 * rank; then for every pi in S_n and every k-subset of positions, contract
 * (drop the k positions and their k values, close the gaps preserving
 * cyclic/linear order on both circles), rank the (n-k)-permutation and look
 * up its precomputed I; track the minimum excess over subsets, then the
 * maximum of that minimum over all pi.
 *
 * Usage: h13i_reduction n k   (n = 4..10, k = 1..3)
 * Version h13i_reduction-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 12

static long factorial_tab[MAXN + 1];

static void init_factorials(void) {
    factorial_tab[0] = 1;
    for (int i = 1; i <= MAXN; i++) factorial_tab[i] = factorial_tab[i - 1] * i;
}

static int inversions(const int *w, int n) {
    int inv = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (w[i] > w[j]) inv++;
    return inv;
}

static int I_of(const int *pi, int n) {
    int best = 1 << 30;
    int w[MAXN];
    for (int q = 0; q < n; q++) {
        for (int c = 0; c < n; c++) {
            int shift = q + 1 - c;
            for (int j = 0; j < n; j++) {
                int v = (pi[(q + 1 + j) % n] - shift) % n;
                if (v < 0) v += n;
                w[j] = v;
            }
            int inv = inversions(w, n);
            if (inv < best) best = inv;
        }
    }
    return best;
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

static long rank_of(const int *a, int m) {
    long r = 0;
    for (int i = 0; i < m; i++) {
        int less = 0;
        for (int j = i + 1; j < m; j++) if (a[j] < a[i]) less++;
        r += (long)less * factorial_tab[m - 1 - i];
    }
    return r;
}

static int f_bound(int m) {
    if (m <= 1) return 0;
    return ((m - 1) * (m - 1)) / 4;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s n k\n", argv[0]); return 2; }
    int n = atoi(argv[1]);
    int k = atoi(argv[2]);
    if (n < 4 || n > MAXN || k < 1 || k > 3 || n - k < 1) return 2;
    init_factorials();

    int m = n - k;
    long total_m = factorial_tab[m];
    unsigned char *Ism = malloc(total_m);
    if (!Ism) { fprintf(stderr, "alloc failed\n"); return 1; }
    {
        int rho[MAXN];
        for (int i = 0; i < m; i++) rho[i] = i;
        long idx = 0;
        do {
            Ism[idx] = (unsigned char)(m >= 2 ? I_of(rho, m) : 0);
            idx++;
        } while (next_perm(rho, m));
        fprintf(stderr, "n=%d k=%d: precomputed I for %ld permutations of size %d\n", n, k, total_m, m);
    }

    long combo[3][20]; /* subsets of size k over n positions, enumerated on the fly */
    int subset[3];
    int pi[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;
    long total_n = factorial_tab[n];
    long rank = 0;
    int worst = -1;
    long worst_rank = -1;
    long counterex = 0;
    int budget = f_bound(n) - f_bound(m);

    do {
        int In = I_of(pi, n);
        int best_excess = 1 << 30;
        /* enumerate all C(n,k) subsets of positions 0..n-1 */
        if (k == 1) {
            for (int a = 0; a < n; a++) {
                subset[0] = a;
                int positions[MAXN], np = 0;
                for (int p = 0; p < n; p++) if (p != subset[0]) positions[np++] = p;
                int dropped[3] = { pi[subset[0]], -1, -1 };
                int contracted[MAXN];
                for (int t = 0; t < np; t++) {
                    int val = pi[positions[t]];
                    int less = 0;
                    for (int d = 0; d < 1; d++) if (dropped[d] < val) less++;
                    contracted[t] = val - less;
                }
                long r2 = rank_of(contracted, m);
                int I2 = Ism[r2];
                int excess = In - I2;
                if (excess < best_excess) best_excess = excess;
            }
        } else if (k == 2) {
            for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) {
                int positions[MAXN], np = 0;
                for (int p = 0; p < n; p++) if (p != a && p != b) positions[np++] = p;
                int dropped[2] = { pi[a], pi[b] };
                if (dropped[0] > dropped[1]) { int t = dropped[0]; dropped[0] = dropped[1]; dropped[1] = t; }
                int contracted[MAXN];
                for (int t = 0; t < np; t++) {
                    int val = pi[positions[t]];
                    int less = (dropped[0] < val) + (dropped[1] < val);
                    contracted[t] = val - less;
                }
                long r2 = rank_of(contracted, m);
                int I2 = Ism[r2];
                int excess = In - I2;
                if (excess < best_excess) best_excess = excess;
            }
        } else { /* k == 3 */
            for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) for (int c2 = b + 1; c2 < n; c2++) {
                int positions[MAXN], np = 0;
                for (int p = 0; p < n; p++) if (p != a && p != b && p != c2) positions[np++] = p;
                int dropped[3] = { pi[a], pi[b], pi[c2] };
                for (int x = 0; x < 3; x++) for (int y = x + 1; y < 3; y++)
                    if (dropped[x] > dropped[y]) { int t = dropped[x]; dropped[x] = dropped[y]; dropped[y] = t; }
                int contracted[MAXN];
                for (int t = 0; t < np; t++) {
                    int val = pi[positions[t]];
                    int less = (dropped[0] < val) + (dropped[1] < val) + (dropped[2] < val);
                    contracted[t] = val - less;
                }
                long r2 = rank_of(contracted, m);
                int I2 = Ism[r2];
                int excess = In - I2;
                if (excess < best_excess) best_excess = excess;
            }
        }
        if (best_excess > worst) { worst = best_excess; worst_rank = rank; }
        if (best_excess > budget) counterex++;
        rank++;
    } while (next_perm(pi, n));

    printf("{\"n\": %d, \"k\": %d, \"m\": %d, \"f_n\": %d, \"f_m\": %d, \"budget\": %d, "
           "\"worst_min_excess\": %d, \"worst_rank\": %ld, \"counterexamples\": %ld, \"total\": %ld}\n",
           n, k, m, f_bound(n), f_bound(m), budget, worst, worst_rank, counterex, total_n);
    free(Ism);
    (void)combo;
    return 0;
}
