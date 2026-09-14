/* line_I_max.c -- exhaustively compute max_pi I(pi) for permutations of n,
 * where I(pi) = min over all n^2 double-cuts (q, c) of inv(w) (line_model.py).
 * Unlike experiments/line_profile.c this does not read/need a distance table,
 * so it is cheap enough to extend the H13-I exhaustive verification range
 * (previously 4 <= n <= 10, docs/notes/h13_line_model.md) to n = 11 (session 9).
 * Usage: line_I_max n
 * Output (stdout): one JSON object with n, total permutations, threshold
 * floor((n-1)^2/4), the observed max_pi I(pi), how many pi attain it, and one
 * example. Version line_I_max-1.0.
 */
#include <stdio.h>
#include <stdlib.h>

#define MAXN 13

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
    int pi[MAXN], w[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;
    long total = 0;
    int maxI = -1;
    int argmax_pi[MAXN];
    long extremal_count = 0;
    int thresh = (n - 1) * (n - 1) / 4;
    do {
        int best = 999;
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
                if (best == 0) goto done_cuts;
            }
        }
        done_cuts:;
        if (best > maxI) {
            maxI = best;
            extremal_count = 1;
            for (int i = 0; i < n; i++) argmax_pi[i] = pi[i];
        } else if (best == maxI) {
            extremal_count++;
        }
        if (best > thresh) {
            fprintf(stderr, "COUNTEREXAMPLE: I=%d exceeds threshold %d at rank %ld\n", best, thresh, total);
        }
        total++;
    } while (next_perm(pi, n));
    printf("{\"n\": %d, \"total\": %ld, \"threshold\": %d, \"max_I\": %d, \"extremal_count\": %ld, \"argmax_pi\": [", n, total, thresh, maxI, extremal_count);
    for (int i = 0; i < n; i++) printf("%d%s", argmax_pi[i], i + 1 < n ? "," : "");
    printf("]}\n");
    return 0;
}
