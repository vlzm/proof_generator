/* h13i_rotation_only.c -- exhaustive check of a restricted one-parameter family
 * of double cuts for H13-I (docs/notes/h13_line_model.md, session 9).
 *
 * Fix the value cut e = 0 (equivalently: fix the position cut q = 0 and let
 * the value cut range, by the identity inv(w) = inv(w^{-1}) applied to
 * pi rotated -- see report), so the only freedom left is the position cut q,
 * i.e. the n cyclic ROTATIONS of the sequence pi(0..n-1) itself (no value
 * relabelling at all):
 *
 *   Q(pi) = min_{q in Z_n} inv(rotate_q(pi))       rotate_q(pi)[j] = pi[(q+j) % n]
 *
 * This is I(pi) restricted to the n cuts {(q, c) : c = q+1} instead of all n^2.
 * inv(rotate_{q+1}(pi)) - inv(rotate_q(pi)) = (n-1) - 2k, where k is the number
 * of remaining entries smaller than the entry leaving the front (the general
 * single-point-flip rule, PLAN/h13 note assignment); this gives an O(n) update
 * per rotation, O(n) per permutation, so all n! permutations are cheap even at
 * n = 11 (39 916 800 permutations).
 *
 * Output: for each n, max_pi Q(pi), the count of maximizers, one maximizer,
 * and the gap against floor((n-1)^2/4).
 *
 * Usage: h13i_rotation_only n
 * Version h13i_rotation_only-1.0.
 */
#include <stdio.h>
#include <stdlib.h>

#define MAXN 16

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
    int w[MAXN], pi[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;
    long worst = -1, worst_count = 0;
    int argworst[MAXN];
    long total = 0;
    do {
        total++;
        for (int i = 0; i < n; i++) w[i] = pi[i];
        int inv = 0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                if (w[i] > w[j]) inv++;
        int best = inv;
        for (int r = 1; r < n; r++) {
            int k = 0;
            for (int j = 1; j < n; j++) if (w[j] < w[0]) k++;
            inv += (n - 1 - 2 * k);
            int w0 = w[0];
            for (int j = 0; j < n - 1; j++) w[j] = w[j + 1];
            w[n - 1] = w0;
            if (inv < best) best = inv;
        }
        if (best > worst) { worst = best; worst_count = 1; for (int i = 0; i < n; i++) argworst[i] = pi[i]; }
        else if (best == worst) worst_count++;
    } while (next_perm(pi, n));
    int thresh = (n - 1) * (n - 1) / 4;
    printf("{\"n\": %d, \"total\": %ld, \"thresh_floor_(n-1)^2/4\": %d, \"max_Q\": %ld, "
           "\"gap\": %ld, \"count_maximizers\": %ld, \"one_argmax\": [",
           n, total, thresh, worst, worst - thresh, worst_count);
    for (int i = 0; i < n; i++) printf("%d%s", argworst[i], i + 1 < n ? ", " : "");
    printf("]}\n");
    return 0;
}
