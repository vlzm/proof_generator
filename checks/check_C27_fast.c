/* check_C27_fast.c -- exhaustive check of H10 / C27 (first part):
 *
 *     min_c [ 2 F_c - S_c + R_c ] <= B_n + 1      for all pi in S_n,
 *
 * where f_c(i) = pi(i) + c mod n, F_c / S_c are the N1 section-1 sums over the
 * nontrivial cycles of f_c and R_c is the length of a shortest walk on Z_n from
 * 0 through every carrier of those cycles to c (carrier-route variant, rule 16).
 *
 * The arithmetic core is checks/carrier_core.h, an independent reimplementation
 * of constructions/strict_upper.py (strict_upper-1.0); checks/check_C27.py
 * cross-checks the two on all permutations for 4 <= n <= 7.
 *
 * Build:  gcc -O2 -o checks/check_C27_fast checks/check_C27_fast.c
 * Run:    ./checks/check_C27_fast <n> [examples_out.txt]
 *
 * Prints B_n, the histogram of min_c[2F_c-S_c+R_c] - B_n over all n!
 * permutations, its maximum, the number of permutations with min_c >= B_n,
 * and the maximum over pi of n*(mean_c value - B_n).  With examples_out.txt
 * every pi with min_c value >= B_n is written there ("excess pi...").
 */

#include <string.h>
#include <time.h>
#include "carrier_core.h"

/* -f mode: read lines "n a0 a1 ... a_{n-1}" from a file and print, for each,
 * the per-shift quadruples "c F_c S_c R_c value" -- used by checks/check_C27.py
 * to cross-check this core against constructions/strict_upper.py.            */
static int run_file(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) { perror("open"); return 1; }
    int m;
    while (fscanf(f, "%d", &m) == 1) {
        int pi[MAXN];
        if (m < 2 || m > MAXN) { fprintf(stderr, "bad n\n"); return 1; }
        for (int i = 0; i < m; i++) if (fscanf(f, "%d", &pi[i]) != 1) { fprintf(stderr, "bad perm\n"); return 1; }
        cc_set_n(m);
        for (int c = 0; c < m; c++) {
            int F, S, R;
            int v = cc_shift_value(pi, c, &F, &S, &R, NULL);
            printf("%d %d %d %d %d%s", c, F, S, R, v, c == m - 1 ? "\n" : " ");
        }
    }
    fclose(f);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s n [examples_out]   |   %s -f perms.txt\n", argv[0], argv[0]); return 1; }
    if (argv[1][0] == '-' && argv[1][1] == 'f') {
        if (argc < 3) { fprintf(stderr, "usage: %s -f perms.txt\n", argv[0]); return 1; }
        return run_file(argv[2]);
    }
    int n = atoi(argv[1]);
    if (n < 2 || n > MAXN) { fprintf(stderr, "n out of range\n"); return 1; }
    cc_set_n(n);
    int Bn = n * (n - 1) / 2;
    FILE *ex = NULL;
    if (argc > 2) { ex = fopen(argv[2], "w"); if (!ex) { perror("open"); return 1; } }

    int pi[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;
    long long count = 0, atleast0 = 0;
    const int hist_lo = -100, hist_hi = 100;
    long long *hist = calloc(hist_hi - hist_lo + 1, sizeof(long long));
    int maxexc = -10000;
    long long maxmean_num = -1;          /* max over pi of sum_c value - n*B_n */
    int maxmean_pi[MAXN];
    clock_t t0 = clock();

    for (;;) {
        long long sum = 0;
        int best = cc_min_over_c(pi, &sum, NULL);
        int exc = best - Bn;
        int he = exc < hist_lo ? hist_lo : (exc > hist_hi ? hist_hi : exc);
        hist[he - hist_lo]++;
        if (exc > maxexc) maxexc = exc;
        if (exc >= 0) {
            atleast0++;
            if (ex) {
                fprintf(ex, "%d", exc);
                for (int i = 0; i < n; i++) fprintf(ex, " %d", pi[i]);
                fprintf(ex, "\n");
            }
        }
        if (sum - (long long)n * Bn > maxmean_num) {
            maxmean_num = sum - (long long)n * Bn;
            memcpy(maxmean_pi, pi, sizeof(int) * n);
        }
        count++;
        int i = n - 2;                            /* next permutation, lex order */
        while (i >= 0 && pi[i] >= pi[i + 1]) i--;
        if (i < 0) break;
        int j = n - 1;
        while (pi[j] <= pi[i]) j--;
        int tmp = pi[i]; pi[i] = pi[j]; pi[j] = tmp;
        for (int a = i + 1, b = n - 1; a < b; a++, b--) { tmp = pi[a]; pi[a] = pi[b]; pi[b] = tmp; }
    }
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("n=%d B_n=%d perms=%lld max(min_c - B_n)=%d count(min_c >= B_n)=%lld time=%.1fs\n",
           n, Bn, count, maxexc, atleast0, secs);
    printf("hist(min_c - B_n):");
    for (int v = hist_lo; v <= hist_hi; v++) if (hist[v - hist_lo]) printf(" %d:%lld", v, hist[v - hist_lo]);
    printf("\n");
    printf("max_pi n*(mean_c - B_n) = %lld  (n=%d) at pi =", maxmean_num, n);
    for (int i = 0; i < n; i++) printf(" %d", maxmean_pi[i]);
    printf("\n");
    if (ex) fclose(ex);
    free(hist);
    return 0;
}
