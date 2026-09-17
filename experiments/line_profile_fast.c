/* line_profile_fast.c -- max over all pi in S_n of I(pi) = min_{q,c} inv(w)
 * (see experiments/line_model.py / experiments/line_profile.c for the
 * definition of the double cut (q, c) and the relabelled line w).
 *
 * Does not need a distance table: I(pi) is a property of pi alone, not of
 * d(pi), so this checks H13-I in isolation and can go further than
 * line_profile.c (which joins against dist_n{n}.bin and is capped by table
 * availability).
 *
 * Speed-up over the brute force (line_profile.c: recompute inversions from
 * scratch for every one of the n^2 cuts, O(n^2) each -> O(n^4) per pi):
 * for a fixed value-shift s, moving the first element x_0 of the line to the
 * end (i.e. advancing the position-cut q by one) changes the inversion count
 * by exactly (n - 1 - 2*x_0), because x_0 was compared as "first" against the
 * other n-1 elements (contributing an inversion for each of the x_0 elements
 * smaller than it) and becomes "last" (contributing an inversion for each of
 * the n-1-x_0 elements larger than it). This turns the sweep over all n
 * position-cuts, for fixed s, into one O(n) prefix-sum walk after a single
 * O(n^2) inversion count, i.e. O(n^2) per shift s and O(n^3) per pi (down
 * from O(n^4)); cross-checked against the brute force (all n^2 cuts,
 * inversions recomputed each time) for 4 <= n <= 9 with --check.
 *
 * Usage: line_profile_fast n [--check]
 * Version line_profile_fast-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXN 14

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

static int I_fast(const int *pi, int n, int *x) {
    int best = 1 << 30;
    for (int s = 0; s < n; s++) {
        for (int i = 0; i < n; i++) {
            int v = pi[i] - s;
            v %= n; if (v < 0) v += n;
            x[i] = v;
        }
        int cur = inversions(x, n);
        if (cur < best) best = cur;
        for (int k = 0; k < n - 1; k++) {
            cur += (n - 1 - 2 * x[k]);
            if (cur < best) best = cur;
        }
    }
    return best;
}

static int I_brute(const int *pi, int n) {
    int best = 1 << 30;
    int w[MAXN];
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
        }
    }
    return best;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s n [--check]\n", argv[0]); return 2; }
    int n = atoi(argv[1]);
    int check = argc >= 3 && !strcmp(argv[2], "--check");
    if (n < 2 || n > MAXN) return 2;

    int pi[MAXN], xscr[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;

    long total = 1;
    for (int k = 2; k <= n; k++) total *= k;

    int bound = ((n - 1) * (n - 1)) / 4;
    int worst = -1;
    long worst_rank = -1;
    long count_at_bound = 0;
    long rank = 0;

    clock_t t0 = clock();
    do {
        int I = I_fast(pi, n, xscr);
        if (check) {
            int Ib = I_brute(pi, n);
            if (Ib != I) {
                printf("MISMATCH at rank %ld: fast=%d brute=%d perm=[", rank, I, Ib);
                for (int i = 0; i < n; i++) printf("%d ", pi[i]);
                printf("]\n");
                return 1;
            }
        }
        if (I > worst) { worst = I; worst_rank = rank; }
        if (I == bound) count_at_bound++;
        rank++;
    } while (next_perm(pi, n));
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;

    printf("{\"n\": %d, \"count\": %ld, \"bound_floor((n-1)^2/4)\": %d, \"max_I\": %d, "
           "\"max_I_minus_bound\": %d, \"argmax_rank\": %ld, \"count_at_bound\": %ld, "
           "\"checked_against_brute\": %s, \"seconds\": %.3f}\n",
           n, total, bound, worst, worst - bound, worst_rank, count_at_bound,
           check ? "true" : "false", secs);
    return 0;
}
