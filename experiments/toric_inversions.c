/* toric_inversions.c — exhaustive study of the toric inversion minimum
 * I(pi) = min over the n^2 double cuts (a, s) of inv(w), w_j = (pi(a+j) - s) mod n
 * (H13-I: I(pi) <= floor((n-1)^2/4)), together with the quantities of the
 * toric identity found in session 9:
 *
 *   inv(a,s) = Q(a,s)/n + c(pi),
 *   Q(a,s)   = sum_j (j - w_j)^2,
 *   c(pi)    = (1/n^2) * sum_{i<k} ( T(k-i) - T(pi(k)-pi(i)) )^2,   T(u) = u mod n.
 *
 * Reported per n (all permutations in lexicographic order):
 *   max_pi I(pi), the number of argmax and whether every argmax is a reflection;
 *   max_pi ( I(pi) - min_r J(r) )  — the wrap-around ("misalignment") obstruction,
 *     J(r) = c + (1/n) sum_i D(d_i - r)^2, d_i = (i-pi(i)) mod n, D = cyclic distance
 *     (J(r) <= inv(a,s) for every cut with (a-s) mod n = r; min_r J(r) <= bound is
 *      proved, see docs/proofs/C37_toric_identity.md);
 *   max_pi sum_a min_s inv(a,s)   — conjecture P' (<= n * bound);
 *   max_pi sum_s min_a inv(a,s)   — the transposed form;
 *   max_pi min_t sum_a inv(a, t-a)      — anti-diagonal average (refuted);
 *   max_pi min_g sum_a inv(a, pi(a)+g)  — shifted data-cut average (refuted).
 *
 * All rational quantities are kept scaled by n^2 (n2c = n^2*c, n2J = n^2*J).
 *
 * Usage: toric_inversions n [n_max]      Version toric_inversions-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 12

static int n;

static int inv_bruteforce(const int *w) {
    int c = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (w[i] > w[j]) c++;
    return c;
}

static int next_perm(int *a, int m) {
    int i = m - 2;
    while (i >= 0 && a[i] >= a[i + 1]) i--;
    if (i < 0) return 0;
    int j = m - 1;
    while (a[j] <= a[i]) j--;
    int t = a[i]; a[i] = a[j]; a[j] = t;
    for (int l = i + 1, r = m - 1; l < r; l++, r--) { t = a[l]; a[l] = a[r]; a[r] = t; }
    return 1;
}

/* inv[a][s] filled by the two increment recursions
 *   inv(a+1,s) - inv(a,s) = (n-1) - 2*((pi(a)-s) mod n)
 *   inv(a,s+1) - inv(a,s) = (n-1) - 2*((pi^{-1}(s)-a) mod n)
 * The recursions are verified against the brute force below for every n. */
static int tab[MAXN][MAXN];

static void fill_table(const int *pi, const int *ipi) {
    int w[MAXN];
    for (int j = 0; j < n; j++) w[j] = pi[j];
    tab[0][0] = inv_bruteforce(w);
    for (int s = 0; s + 1 < n; s++) {
        int lam = ipi[s] % n;
        tab[0][s + 1] = tab[0][s] + (n - 1) - 2 * lam;
    }
    for (int s = 0; s < n; s++)
        for (int a = 0; a + 1 < n; a++) {
            int rho = (pi[a] - s) % n; if (rho < 0) rho += n;
            tab[a + 1][s] = tab[a][s] + (n - 1) - 2 * rho;
        }
}

static void verify_table(const int *pi) {
    int w[MAXN];
    for (int a = 0; a < n; a++)
        for (int s = 0; s < n; s++) {
            for (int j = 0; j < n; j++) {
                int v = (pi[(a + j) % n] - s) % n; if (v < 0) v += n;
                w[j] = v;
            }
            if (inv_bruteforce(w) != tab[a][s]) {
                fprintf(stderr, "recursion mismatch n=%d a=%d s=%d\n", n, a, s);
                exit(1);
            }
        }
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s n [n_max]\n", argv[0]); return 2; }
    int n0 = atoi(argv[1]);
    int n1 = (argc >= 3) ? atoi(argv[2]) : n0;
    printf("{\"tool\": \"toric_inversions-1.0\", \"runs\": [\n");
    for (n = n0; n <= n1; n++) {
        if (n < 3 || n > MAXN) return 2;
        int bound = ((n - 1) * (n - 1)) / 4;
        long total = 1;
        for (int k = 2; k <= n; k++) total *= k;

        int pi[MAXN], ipi[MAXN];
        for (int i = 0; i < n; i++) pi[i] = i;

        int maxI = -1; long argI = -1, cntI = 0, cntI_refl = 0, violations = 0;
        long long maxGapJ = -1; long argGapJ = -1;      /* n^2 * (I - min_r J) */
        long maxRowSum = -1, argRow = -1;               /* sum_a min_s inv */
        long maxColSum = -1;
        long maxAnti = -1, argAnti = -1;                /* min_t sum_a inv(a,t-a) */
        long maxShift = -1, argShift = -1;              /* min_g sum_a inv(a,pi(a)+g) */
        long rank = 0;
        long verified_cnt = 0;

        do {
            for (int i = 0; i < n; i++) ipi[pi[i]] = i;
            fill_table(pi, ipi);
            /* independent check of the two increment recursions against the
             * definition: on the first 100 permutations and then every 997-th */
            if (rank < 100 || rank % 997 == 0) { verify_table(pi); verified_cnt++; }

            int I = tab[0][0];
            for (int a = 0; a < n; a++)
                for (int s = 0; s < n; s++)
                    if (tab[a][s] < I) I = tab[a][s];
            if (I > bound) violations++;
            if (I > maxI) { maxI = I; argI = rank; cntI = 0; cntI_refl = 0; }
            if (I == maxI) {
                cntI++;
                /* reflection test: pi(i) + i constant mod n */
                int h = (pi[0] + 0) % n, isrefl = 1;
                for (int i = 1; i < n; i++) if ((pi[i] + i) % n != h) { isrefl = 0; break; }
                if (isrefl) cntI_refl++;
            }

            /* c scaled: n2c = n^2 * c = sum over unordered pairs of (T(dk)-T(dv))^2 */
            long long n2c = 0;
            for (int i = 0; i < n; i++)
                for (int k = i + 1; k < n; k++) {
                    int p = (k - i) % n;
                    int q = (pi[k] - pi[i]) % n; if (q < 0) q += n;
                    n2c += (long long)(p - q) * (p - q);
                }
            /* min_r J(r):  n^2 J(r) = n2c + n * sum_i D(d_i - r)^2 */
            long long minJ = -1;
            for (int r = 0; r < n; r++) {
                long long sq = 0;
                for (int i = 0; i < n; i++) {
                    int d = (i - pi[i] - r) % n; if (d < 0) d += n;
                    int D = d < n - d ? d : n - d;
                    sq += (long long)D * D;
                }
                long long J = n2c + (long long)n * sq;
                if (minJ < 0 || J < minJ) minJ = J;
            }
            long long gap = (long long)n * n * I - minJ;
            if (gap > maxGapJ) { maxGapJ = gap; argGapJ = rank; }

            long rowsum = 0, colsum = 0;
            for (int a = 0; a < n; a++) {
                int m = tab[a][0];
                for (int s = 1; s < n; s++) if (tab[a][s] < m) m = tab[a][s];
                rowsum += m;
            }
            for (int s = 0; s < n; s++) {
                int m = tab[0][s];
                for (int a = 1; a < n; a++) if (tab[a][s] < m) m = tab[a][s];
                colsum += m;
            }
            if (rowsum > maxRowSum) { maxRowSum = rowsum; argRow = rank; }
            if (colsum > maxColSum) maxColSum = colsum;

            long bestAnti = -1, bestShift = -1;
            for (int t = 0; t < n; t++) {
                long sa = 0, sb = 0;
                for (int a = 0; a < n; a++) {
                    sa += tab[a][((t - a) % n + n) % n];
                    sb += tab[a][(pi[a] + t) % n];
                }
                if (bestAnti < 0 || sa < bestAnti) bestAnti = sa;
                if (bestShift < 0 || sb < bestShift) bestShift = sb;
            }
            if (bestAnti > maxAnti) { maxAnti = bestAnti; argAnti = rank; }
            if (bestShift > maxShift) { maxShift = bestShift; argShift = rank; }

            rank++;
        } while (next_perm(pi, n));

        printf("  {\"n\": %d, \"perms\": %ld, \"bound\": %d, \"max_I\": %d, "
               "\"H13I_violations\": %ld, \"argmax_I_rank\": %ld, \"argmax_count\": %ld, "
               "\"argmax_reflections\": %ld,\n"
               "   \"max_I_minus_minJ_x_n2\": %lld, \"n2\": %d, \"argmax_gapJ_rank\": %ld,\n"
               "   \"max_sum_a_min_s_inv\": %ld, \"n_times_bound\": %d, \"argmax_rowsum_rank\": %ld,\n"
               "   \"max_sum_s_min_a_inv\": %ld,\n"
               "   \"max_min_t_antidiag\": %ld, \"argmax_antidiag_rank\": %ld,\n"
               "   \"max_min_g_shifted_data\": %ld, \"argmax_shifted_rank\": %ld,\n"
               "   \"recursion_verified_perms\": %ld}%s\n",
               n, total, bound, maxI, violations, argI, cntI, cntI_refl,
               maxGapJ, n * n, argGapJ,
               maxRowSum, n * bound, argRow,
               maxColSum,
               maxAnti, argAnti,
               maxShift, argShift, verified_cnt,
               n == n1 ? "" : ",");
        fflush(stdout);
    }
    printf("]}\n");
    return 0;
}
