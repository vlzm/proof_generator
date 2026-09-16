/* line_toric_scan.c — exhaustive scan of H13-I over toric classes.
 *
 * I(pi) = min over the n^2 double cuts (position cut q, value cut c) of the
 * inversion count of the line w_j = (pi(q+1+j) - (q+1-c)) mod n.  I is constant
 * on a toric class (rotations of positions and of values), and every toric class
 * contains a line with w_0 = 0, so it is enough to enumerate the (n-1)!
 * permutations w with w_0 = 0.
 *
 * For each such representative the whole n x n grid of inversion counts is built
 * with O(1) updates:
 *   inv(a, b+1) = inv(a, b) + (n - 1 - 2p),  p = position of value b in line(a,b);
 *   inv(a+1, 0) = inv(a, 0) + (n - 1 - 2*w_a).
 *
 * The same sweep also tests the rule of CLAIMS C39 (H13-S): among the cuts
 * maximising sum_j j*w_j (equivalently minimising the Spearman distance
 * sum_j (w_j - j)^2) take the worst inversion count; the weak form asks that it
 * be at most floor((n-1)^2/4), the strong form that it equal I(pi).  The score
 * is updated in O(1) as well: sum_j j*w_j changes by n*p - n(n-1)/2 (value cut)
 * and by n*w_a - n(n-1)/2 (position cut).
 *
 * Output (stdout, one JSON object): n, number of classes scanned, max I,
 * floor((n-1)^2/4), the number of representatives attaining the max, the first
 * few of them, the histogram of I, and the two H13-S failure counters.
 *
 * Usage: line_toric_scan n
 * Version line_toric_scan-1.1 (1.0 + the H13-S counters).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
    if (argc < 2) { fprintf(stderr, "usage: %s n\n", argv[0]); return 2; }
    int n = atoi(argv[1]);
    if (n < 3 || n > MAXN) return 2;

    int w[MAXN], pos[MAXN];
    long hist[MAXI];
    memset(hist, 0, sizeof(hist));
    long classes = 0, nmax = 0, hs_weak_fail = 0, hs_strong_fail = 0;
    int best = -1;
    long shown = 0;
    char examples[4096]; examples[0] = 0;

    for (int i = 0; i < n; i++) w[i] = i;   /* w[0] = 0 fixed, permute w[1..n-1] */

    do {
        classes++;
        for (int i = 0; i < n; i++) pos[w[i]] = i;
        /* inv(0,0) */
        int inv0 = 0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                if (w[i] > w[j]) inv0++;
        long s0 = 0;
        for (int j = 0; j < n; j++) s0 += (long)j * w[j];
        const long T = (long)n * (n - 1) / 2;
        int cur = inv0, mn = inv0;
        long curS = s0, bestS = s0;
        int bestSI = inv0;
        for (int a = 0; a < n; a++) {
            int v = cur; long vS = curS;
            for (int b = 0; b < n; b++) {
                if (vS > bestS) { bestS = vS; bestSI = v; }
                else if (vS == bestS && v > bestSI) bestSI = v;
                int p = pos[b] - a; if (p < 0) p += n;
                v += n - 1 - 2 * p;
                vS += (long)n * p - T;
                if (v < mn) mn = v;
            }
            /* v is now inv(a,0) again (full cycle of value shifts) */
            cur += n - 1 - 2 * w[a];
            curS += (long)n * w[a] - T;
        }
        if (bestSI > (n - 1) * (n - 1) / 4) hs_weak_fail++;
        if (bestSI != mn) hs_strong_fail++;
        hist[mn]++;
        if (mn > best) {
            best = mn; nmax = 1; shown = 0; examples[0] = 0;
        } else if (mn == best) nmax++;
        if (mn == best && shown < 5) {
            char buf[128]; int k = 0;
            k += snprintf(buf + k, sizeof(buf) - k, "%s[", shown ? ", " : "");
            for (int i = 0; i < n; i++) k += snprintf(buf + k, sizeof(buf) - k, "%s%d", i ? "," : "", w[i]);
            k += snprintf(buf + k, sizeof(buf) - k, "]");
            if (strlen(examples) + k < sizeof(examples) - 1) strcat(examples, buf);
            shown++;
        }
    } while (next_perm(w + 1, n - 1));

    printf("{\"version\": \"line_toric_scan-1.1\", \"n\": %d, \"classes_scanned\": %ld,\n", n, classes);
    printf(" \"max_I\": %d, \"floor((n-1)^2/4)\": %d, \"argmax_count\": %ld,\n", best, (n - 1) * (n - 1) / 4, nmax);
    printf(" \"argmax_examples\": [%s],\n", examples);
    printf(" \"H13S_weak_failures\": %ld, \"H13S_strong_failures\": %ld,\n",
           hs_weak_fail, hs_strong_fail);
    printf(" \"hist_I\": {");
    int first = 1;
    for (int i = 0; i < MAXI; i++)
        if (hist[i]) { printf("%s\"%d\": %ld", first ? "" : ", ", i, hist[i]); first = 0; }
    printf("}}\n");
    return 0;
}
