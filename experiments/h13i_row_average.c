/* h13i_row_average.c — session 9 probe for H13-I (docs/notes/h13_line_model.md
 * S7): exhaustive check of the new conjecture
 *
 *   sum_{q=0}^{n-1} I(q) <= n * floor((n-1)^2/4),
 *
 * where I(q) = min_s inv(q,s) is the row-minimum of the double-cut model
 * (position cut q, value cut s) already used in experiments/line_model.py
 * and experiments/line_profile.c: for the row a_j = pi((q+1+j) mod n),
 * j = 0..n-1, and w_j = (a_j - s) mod n, inv(q,s) = #{j<k : w_j > w_k}.
 *
 * This is strictly stronger than H13-I itself (I(pi) = min_q I(q) <=
 * average_q I(q), so the inequality above implies I(pi) <=
 * floor((n-1)^2/4) by pigeonhole) and is NOT the double average over all
 * n^2 cuts that is already known to fail (identity averages inv(q,s) over
 * s for FIXED q first, i.e. it is an average of true per-row MINIMA, not
 * an average of inv itself).
 *
 * Closed form used for I(q) (derived and checked against brute force over
 * all (q,s) in experiments/h13i_probe.py, part 2):
 *   rho = a^{-1} (position, within the row, of each value 0..n-1)
 *   inv_a = standard inversions of the row a
 *   R(s) = rho(0) + rho(1) + ... + rho(s-1)   (s = 0..n-1, R(0) = 0)
 *   M(s) = 2*R(s) - s*(n-1)
 *   I(q) = inv_a - max_{0<=s<=n-1} M(s)
 * (M(s) is (twice) the excess, over the s smallest values, of their actual
 * position-sum in the row over the "balanced" position-sum s*(n-1)/2; see
 * docs/notes/h13_line_model.md S7.2 for the derivation from the recursion
 * inv(q,s+1) - inv(q,s) = (n-1) - 2*((pi^{-1}(s) - q - 1) mod n).)
 *
 * Usage: h13i_row_average n
 * Output: one JSON line with n, the target bound, n*bound, the number of
 * permutations checked (n! exhaustively), the number of violations, the
 * largest excess found (sum_q I(q) - n*bound) and its argmax permutation.
 * Version h13i_row_average-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 14

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
    if (argc < 2) { fprintf(stderr, "usage: %s n\n", argv[0]); return 1; }
    int n = atoi(argv[1]);
    if (n < 2 || n > MAXN) { fprintf(stderr, "n out of range\n"); return 1; }
    int pi[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;

    int bound = ((n - 1) * (n - 1)) / 4;
    long long nbound = (long long)n * bound;

    long long checked = 0, violations = 0;
    long long worst_excess = -1000000000LL;
    int worst_pi[MAXN];
    int row[MAXN], rho[MAXN];

    do {
        long long total = 0;
        for (int q = 0; q < n; q++) {
            for (int j = 0; j < n; j++) row[j] = pi[(q + 1 + j) % n];
            for (int j = 0; j < n; j++) rho[row[j]] = j;
            int inv_a = 0;
            for (int j = 0; j < n; j++)
                for (int k = j + 1; k < n; k++)
                    if (row[j] > row[k]) inv_a++;
            int R = 0, bestM = 0; /* s = 0 gives M(0) = 0 */
            for (int s = 1; s < n; s++) {
                R += rho[s - 1];
                int M = 2 * R - s * (n - 1);
                if (M > bestM) bestM = M;
            }
            total += inv_a - bestM;
        }
        checked++;
        long long excess = total - nbound;
        if (excess > 0) violations++;
        if (excess > worst_excess) {
            worst_excess = excess;
            memcpy(worst_pi, pi, sizeof(int) * n);
        }
    } while (next_perm(pi, n));

    printf("{\"n\": %d, \"bound\": %d, \"n_times_bound\": %lld, \"checked\": %lld, "
           "\"violations\": %lld, \"worst_excess\": %lld, \"worst_pi\": [",
           n, bound, nbound, checked, violations, worst_excess);
    for (int i = 0; i < n; i++) printf("%d%s", worst_pi[i], i + 1 < n ? "," : "");
    printf("]}\n");
    return 0;
}
