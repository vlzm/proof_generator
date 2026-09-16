/* toric_cut_diag_test.c — tests candidate restrictions of the (a,b) search
 * for H13-I (see docs/notes/h13_line_model.md §6): instead of the full n^2
 * cuts, only check cuts (a,b) with some cheap-to-describe relation, and see
 * whether the restricted minimum still meets floor((n-1)^2/4) on every
 * permutation. Uses the O(n^2) recursion from toric_cut_fast.c (D1).
 *
 * Modes (each: exhaustive over all permutations of n, report any pi where
 * the restricted minimum exceeds floor((n-1)^2/4), with a counterexample):
 *   diag n     — b = pi[a] (cut lands exactly on a point of pi), n candidates
 *   diag2 n    — b = pi[a] or b = pi[a]+1 or b=pi[a]-1 (3n candidates)
 *   fixed_a n  — a = 0 only (all b), 1 candidate a, n candidates b (Lemma C style)
 * Usage: toric_cut_diag_test mode n
 * Version toric_cut_diag_test-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 14

static void full_table(const int *pi, const int *inv_pi, int n, int T[MAXN][MAXN]) {
    int inv00 = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (pi[i] > pi[j]) inv00++;
    int Aa[MAXN];
    Aa[0] = inv00;
    for (int a = 0; a < n - 1; a++) {
        int val = pi[a] % n;
        Aa[a + 1] = Aa[a] + (n - 1 - 2 * val);
    }
    for (int a = 0; a < n; a++) {
        T[a][0] = Aa[a];
        int cur = Aa[a];
        for (int b = 0; b < n - 1; b++) {
            int istar = (inv_pi[b] - a) % n;
            if (istar < 0) istar += n;
            cur = cur + (n - 1 - 2 * istar);
            T[a][b + 1] = cur;
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
    if (argc < 3) { fprintf(stderr, "usage: %s mode n\n", argv[0]); return 2; }
    const char *mode = argv[1];
    int n = atoi(argv[2]);
    if (n < 2 || n > MAXN) return 2;
    int bound = ((n - 1) * (n - 1)) / 4;

    int pi[MAXN], inv_pi[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;

    long total = 0, violations = 0;
    int T[MAXN][MAXN];
    do {
        for (int i = 0; i < n; i++) inv_pi[pi[i]] = i;
        full_table(pi, inv_pi, n, T);
        int best = 1000000;
        if (strcmp(mode, "diag") == 0) {
            for (int a = 0; a < n; a++) {
                int b = pi[a];
                if (T[a][b] < best) best = T[a][b];
            }
        } else if (strcmp(mode, "diag2") == 0) {
            for (int a = 0; a < n; a++) {
                for (int d = -1; d <= 1; d++) {
                    int b = ((pi[a] + d) % n + n) % n;
                    if (T[a][b] < best) best = T[a][b];
                }
            }
        } else if (strcmp(mode, "fixed_a") == 0) {
            for (int b = 0; b < n; b++) if (T[0][b] < best) best = T[0][b];
        } else {
            fprintf(stderr, "unknown mode\n"); return 2;
        }
        total++;
        if (best > bound) {
            violations++;
            if (violations <= 5) {
                printf("VIOLATION pi=(");
                for (int i = 0; i < n; i++) printf("%d%s", pi[i], i + 1 < n ? "," : "");
                printf(") restricted_best=%d bound=%d\n", best, bound);
            }
        }
    } while (next_perm(pi, n));
    printf("{\"mode\": \"%s\", \"n\": %d, \"total\": %ld, \"violations\": %ld, \"bound\": %d}\n",
           mode, n, total, violations, bound);
    return 0;
}
