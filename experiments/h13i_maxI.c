/* h13i_maxI.c — exhaustive check of H13-I: for all permutations pi of n,
 * compute I(pi) = min over cut (q,c) of inv(w_{q,c}) (see experiments/line_model.py,
 * experiments/line_profile.c for the exact same definition) and report max_pi I(pi),
 * to compare against floor((n-1)^2/4).  Self-contained: no distance table needed
 * (unlike line_profile.c, which also compares against d(pi) from data/tables/).
 * Used in session 9 to extend the H13-I exhaustive check to n = 11 if the
 * estimated runtime is acceptable (AGENTS.md rule 10).
 *
 * Usage: h13i_maxI n [argmax_only]
 * Version h13i_maxI-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAXN 12

static int inversions(const int *w, int n) {
    int inv = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (w[i] > w[j]) inv++;
    return inv;
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
    if (argc < 2) { fprintf(stderr, "usage: %s n\n", argv[0]); return 2; }
    int n = atoi(argv[1]);
    if (n < 2 || n > MAXN) return 2;
    long total = 1;
    for (int k = 2; k <= n; k++) total *= k;

    int pi[MAXN], w[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;
    long rank = 0;
    int maxI = -1;
    int argmax[MAXN];
    long bound = (long)(n - 1) * (n - 1) / 4;
    long countI[400];
    for (int i = 0; i < 400; i++) countI[i] = 0;
    time_t t0 = time(NULL);
    do {
        int best = 999999;
        for (int q = 0; q < n; q++) {
            int shift0 = q + 1;
            for (int c = 0; c < n; c++) {
                int shift = shift0 - c;
                if (shift < 0) shift += n; /* keep small, mod applied below anyway */
                for (int j = 0; j < n; j++) {
                    int v = (pi[(q + 1 + j) % n] - shift) % n;
                    if (v < 0) v += n;
                    w[j] = v;
                }
                int inv = inversions(w, n);
                if (inv < best) best = inv;
            }
        }
        countI[best]++;
        if (best > maxI) { maxI = best; for (int i = 0; i < n; i++) argmax[i] = pi[i]; }
        rank++;
        if (rank % 20000000 == 0) {
            fprintf(stderr, "  ... %ld / %ld (%.0f s)\n", rank, total, difftime(time(NULL), t0));
        }
    } while (next_perm(pi, n));

    double secs = difftime(time(NULL), t0);
    printf("n=%d total=%ld bound=floor((n-1)^2/4)=%ld max_I=%d %s seconds=%.0f argmax=(", n, total, bound, maxI,
           (long)maxI <= bound ? "OK" : "VIOLATION", secs);
    for (int i = 0; i < n; i++) printf("%d%s", argmax[i], i + 1 < n ? "," : "");
    printf(")\n");
    printf("hist_I:");
    for (int i = 0; i <= maxI; i++) if (countI[i]) printf(" %d:%ld", i, countI[i]);
    printf("\n");
    return (long)maxI <= bound ? 0 : 1;
}
