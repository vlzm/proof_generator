/* h13i_exhaustive.c — exhaustive scan of I(pi) for H13-I, independent
 * re-implementation of I(pi) (see h13i_search.c for the O(n^2) recurrences;
 * this file duplicates them so a bug in one is not shared with the other,
 * per AGENTS.md rule 3). Enumerates all n! permutations in lexicographic
 * order and reports the maximum I(pi) found, compared to floor((n-1)^2/4).
 * Session 8 (line_profile.c) already covers 4<=n<=10 against the distance
 * tables; this tool is used here to extend the same exhaustive check to
 * n=11 with a measured time budget (AGENTS.md rule 10) without needing the
 * distance table (which is not required for I(pi) alone).
 * Usage: h13i_exhaustive n [report_every_seconds]
 * Version h13i_exhaustive-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static int n;
static int pi_arr[16], piinv[16], pos_a[16];

static long inv0_full(const int *row) {
    long c = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (row[i] > row[j]) c++;
    return c;
}

static long I_of_pi(void) {
    for (int v = 0; v < n; v++) piinv[pi_arr[v]] = v;
    long inv_a0 = inv0_full(pi_arr);
    long best = inv_a0, inv_a = inv_a0;
    for (int a = 0; a < n; a++) {
        for (int v = 0; v < n; v++) { int p = piinv[v] - a; if (p < 0) p += n; pos_a[v] = p; }
        long inv_b = inv_a;
        if (inv_b < best) best = inv_b;
        for (int b = 0; b < n - 1; b++) {
            inv_b += (long)(n - 1) - 2 * pos_a[b];
            if (inv_b < best) best = inv_b;
        }
        inv_a += (long)(n - 1) - 2 * pi_arr[a];
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

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s n [report_every_seconds]\n", argv[0]); return 2; }
    n = atoi(argv[1]);
    if (n < 2 || n > 14) { fprintf(stderr, "n out of range\n"); return 2; }
    int report_every = argc >= 3 ? atoi(argv[2]) : 30;
    long bound = ((long)(n - 1) * (n - 1)) / 4;
    for (int i = 0; i < n; i++) pi_arr[i] = i;
    long count = 0;
    long best = -1;
    int best_pi[16];
    long exceed_count = 0;
    time_t start = time(NULL), last_report = start;
    do {
        long v = I_of_pi();
        count++;
        if (v > best) { best = v; for (int i = 0; i < n; i++) best_pi[i] = pi_arr[i]; }
        if (v > bound) exceed_count++;
        time_t now = time(NULL);
        if (now - last_report >= report_every) {
            fprintf(stderr, "n=%d count=%ld best=%ld bound=%ld elapsed=%lds\n", n, count, best, bound, (long)(now - start));
            last_report = now;
        }
    } while (next_perm(pi_arr, n));
    time_t end = time(NULL);
    printf("{\"n\": %d, \"count\": %ld, \"bound\": %ld, \"max_I\": %ld, \"exceed_count\": %ld, \"seconds\": %ld, \"argmax_pi\": [",
           n, count, bound, best, exceed_count, (long)(end - start));
    for (int i = 0; i < n; i++) printf("%d%s", best_pi[i], i + 1 < n ? "," : "");
    printf("]}\n");
    return 0;
}
