/* toric_min_inversions.c — exhaustive maximum of I(pi) over all toric classes.
 *
 * I(pi) = min over the n^2 double cuts (a, s) of the number of inversions of the
 * line w_j = (pi[(a + j) mod n] - s) mod n, j = 0..n-1 (see
 * docs/notes/h13_line_model.md §0; there a = q + 1, s = q + 1 - c).
 * I is constant on the toric class {pi(. + a) - s}, and every class contains a
 * representative with pi[0] = 0 (take a = 0, s = pi(0)), so the scan runs over
 * the (n-1)! representatives with pi[0] = 0 instead of all n! permutations.
 *
 * The whole n x n table inv(a, s) is built with O(1) per cut from
 *   inv(a+1, s) - inv(a, s) = n - 1 - 2 * w_0,   w_0 = (pi[a] - s) mod n,
 *   inv(a, s+1) - inv(a, s) = n - 1 - 2 * p_0,   p_0 = (pi^{-1}(s) - a) mod n
 * (docs/notes/h13_line_model.md, Lemma C), so one permutation costs O(n^2)
 * instead of the O(n^4) of experiments/line_profile.c (which also needs the
 * distance table; this program does not).
 *
 * Output (stdout, one JSON object): n, number of representatives, target
 * floor((n-1)^2/4), max I with the number of representatives attaining it and
 * the first few of them, the histogram of I, the minimum over pi of the number
 * of cuts with inv <= target ("good cuts") with an example, and timing.
 *
 * Usage: toric_min_inversions n [--progress]
 * Version toric_min_inv-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXN 14
#define MAXI 128

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
    if (argc < 2) { fprintf(stderr, "usage: %s n [--progress]\n", argv[0]); return 2; }
    int n = atoi(argv[1]);
    if (n < 3 || n > MAXN) return 2;
    int progress = (argc > 2 && strcmp(argv[2], "--progress") == 0);
    int target = ((n - 1) * (n - 1)) / 4;

    long hist[MAXI];
    memset(hist, 0, sizeof hist);
    int pi[MAXN], pinv[MAXN], row[MAXN];
    int best_pi[8][MAXN];
    int best_pi_cnt = 0;
    int maxI = -1;
    long maxI_count = 0;
    int min_good = 1 << 30;
    int min_good_pi[MAXN];
    long total = 0;
    clock_t t0 = clock();

    for (int i = 0; i < n; i++) pi[i] = i;
    do {
        /* pi[0] = 0 always: next_perm is applied to pi + 1 below. */
        for (int i = 0; i < n; i++) pinv[pi[i]] = i;
        /* inv(0,0) */
        int base = 0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                if (pi[i] > pi[j]) base++;
        int I = 1 << 30, good = 0;
        int cur_a0 = base;
        for (int a = 0; a < n; a++) {
            if (a > 0) cur_a0 += n - 1 - 2 * pi[a - 1];
            int cur = cur_a0;
            row[0] = cur;
            for (int s = 1; s < n; s++) {
                int p0 = pinv[s - 1] - a;
                if (p0 < 0) p0 += n;
                cur += n - 1 - 2 * p0;
                row[s] = cur;
            }
            for (int s = 0; s < n; s++) {
                if (row[s] < I) I = row[s];
                if (row[s] <= target) good++;
            }
        }
        hist[I]++;
        total++;
        if (I > maxI) {
            maxI = I; maxI_count = 1; best_pi_cnt = 0;
            memcpy(best_pi[best_pi_cnt++], pi, sizeof(int) * n);
        } else if (I == maxI) {
            maxI_count++;
            if (best_pi_cnt < 8) memcpy(best_pi[best_pi_cnt++], pi, sizeof(int) * n);
        }
        if (good < min_good) { min_good = good; memcpy(min_good_pi, pi, sizeof(int) * n); }
        if (progress && (total % 10000000L) == 0)
            fprintf(stderr, "  %ld reps, %.0f s\n", total, (double)(clock() - t0) / CLOCKS_PER_SEC);
    } while (next_perm(pi + 1, n - 1));

    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("{\"version\": \"toric_min_inv-1.0\", \"n\": %d, \"reps\": %ld, \"target\": %d,\n",
           n, total, target);
    printf(" \"max_I\": %d, \"max_I_count\": %ld, \"max_I_examples\": [", maxI, maxI_count);
    for (int k = 0; k < best_pi_cnt; k++) {
        printf("%s[", k ? ", " : "");
        for (int i = 0; i < n; i++) printf("%s%d", i ? ", " : "", best_pi[k][i]);
        printf("]");
    }
    printf("],\n \"min_good_cuts\": %d, \"min_good_cuts_pi\": [", min_good);
    for (int i = 0; i < n; i++) printf("%s%d", i ? ", " : "", min_good_pi[i]);
    printf("],\n \"hist_I\": [");
    int first = 1;
    for (int i = 0; i < MAXI; i++) if (hist[i]) {
        printf("%s[%d, %ld]", first ? "" : ", ", i, hist[i]);
        first = 0;
    }
    printf("],\n \"seconds\": %.1f}\n", secs);
    return 0;
}
