/* torus_cut_scan.c — exhaustive scan of the double-cut (torus) model of H13-I.
 *
 * For a permutation pi of Z_n and a double cut (a, b) let
 *     w_j = (pi(a + j) - b) mod n,   j = 0..n-1,
 *     S(a, b) = B_n - 2 * inv(w),    B_n = n(n-1)/2.
 * Then I(pi) = min_{a,b} inv(w) and H13-I says max_{a,b} S(a,b) >= floor(n/2).
 *
 * The scan reports, for every n in [nmin, nmax]:
 *   (1) min over pi of max_{a,b} S(a,b)          (H13-I, equality set counted);
 *   (2) max over pi of sum_a min_b inv(a, b)     (claim A: <= n*floor((n-1)^2/4));
 *   (3) max over lines y of min_b inv(y - b)     (excess over floor((n-1)^2/4):
 *       0,0,0,0,1,1,1,3,3,3,5 for n = 3..13, i.e. one cut alone is not enough);
 *   (4) max over pi and a of (min_b inv(a, b) - floor((n-1)^2/4)), the per-row
 *       excess that claim A has to compensate.
 *
 * Mode "lines" only runs item (3), restricted to y with y_0 = 0 (min_b inv is
 * invariant under a cyclic shift of the values, so this loses nothing).
 *
 * Usage:  torus_cut_scan perms nmin nmax
 *         torus_cut_scan lines nmin nmax
 * Version torus_cut_scan-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXN 13

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

static int inv_line(const int *w, int n) {
    int s = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (w[i] > w[j]) s++;
    return s;
}

/* inv[a][b] for all n^2 cuts, built from inv(0,0) by the exact increments
 *   inv(a+1,b) - inv(a,b) = (n-1) - 2*g,  g = (pi(a) - b) mod n
 *   inv(a,b+1) - inv(a,b) = (n-1) - 2*h,  h = (pi^{-1}(b) - a) mod n     */
static void cut_table(const int *pi, int n, int *inv) {
    int pinv[MAXN], w[MAXN];
    for (int i = 0; i < n; i++) pinv[pi[i]] = i;
    for (int j = 0; j < n; j++) w[j] = pi[j];
    inv[0] = inv_line(w, n);
    for (int b = 0; b + 1 < n; b++) {
        int h = pinv[b];                      /* a = 0 */
        inv[b + 1] = inv[b] + (n - 1) - 2 * h;
    }
    for (int a = 0; a + 1 < n; a++)
        for (int b = 0; b < n; b++) {
            int g = pi[a] - b; if (g < 0) g += n;
            inv[(a + 1) * n + b] = inv[a * n + b] + (n - 1) - 2 * g;
        }
}

static void print_perm(const int *p, int n) {
    printf("[");
    for (int i = 0; i < n; i++) printf(i ? ",%d" : "%d", p[i]);
    printf("]");
}

static void run_perms(int nmin, int nmax) {
    for (int n = nmin; n <= nmax; n++) {
        int B = n * (n - 1) / 2;
        int bound = ((n - 1) * (n - 1)) / 4;
        int half = n / 2;
        int *inv = malloc(sizeof(int) * n * n);
        int p[MAXN], bestS_arg[MAXN], claimA_arg[MAXN], line_arg[MAXN];
        for (int i = 0; i < n; i++) p[i] = i;
        int minmaxS = 1 << 30, minmaxS_cnt = 0;
        long claimA_max = -1; int claimA_cnt = 0;
        int line_max = -1, line_cnt = 0;
        int rowexc_max = -1000;            /* max over pi,a of (min_b inv - bound) */
        clock_t t0 = clock();
        do {
            cut_table(p, n, inv);
            int maxS = -1 << 30;
            long sum_rowmin = 0;
            for (int a = 0; a < n; a++) {
                int rowmin = 1 << 30;
                for (int b = 0; b < n; b++) {
                    int v = inv[a * n + b];
                    if (v < rowmin) rowmin = v;
                }
                sum_rowmin += rowmin;
                int s = B - 2 * rowmin;
                if (s > maxS) maxS = s;
                if (rowmin - bound > rowexc_max) rowexc_max = rowmin - bound;
                if (a == 0 && rowmin > line_max) { line_max = rowmin; line_cnt = 1;
                    memcpy(line_arg, p, sizeof(int) * n); }
                else if (a == 0 && rowmin == line_max) line_cnt++;
            }
            if (maxS < minmaxS) { minmaxS = maxS; minmaxS_cnt = 1;
                memcpy(bestS_arg, p, sizeof(int) * n); }
            else if (maxS == minmaxS) minmaxS_cnt++;
            if (sum_rowmin > claimA_max) { claimA_max = sum_rowmin; claimA_cnt = 1;
                memcpy(claimA_arg, p, sizeof(int) * n); }
            else if (sum_rowmin == claimA_max) claimA_cnt++;
        } while (next_perm(p, n));
        double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
        printf("{\"n\":%d,\"B\":%d,\"bound\":%d,\"half\":%d,", n, B, bound, half);
        printf("\"minmax_S\":%d,\"minmax_S_count\":%d,\"minmax_S_arg\":", minmaxS, minmaxS_cnt);
        print_perm(bestS_arg, n);
        printf(",\"H13I\":%s,", minmaxS >= half ? "true" : "false");
        printf("\"claimA_max\":%ld,\"claimA_target\":%ld,\"claimA_count\":%d,\"claimA_arg\":",
               claimA_max, (long)n * bound, claimA_cnt);
        print_perm(claimA_arg, n);
        printf(",\"claimA\":%s,", claimA_max <= (long)n * bound ? "true" : "false");
        printf("\"line_max\":%d,\"line_excess\":%d,\"line_count\":%d,\"line_arg\":",
               line_max, line_max - bound, line_cnt);
        print_perm(line_arg, n);
        printf(",\"row_excess_max\":%d,\"seconds\":%.1f}\n", rowexc_max, secs);
        fflush(stdout);
        free(inv);
    }
}

/* item (3) alone, over all lines y with y_0 = 0 */
static void run_lines(int nmin, int nmax) {
    for (int n = nmin; n <= nmax; n++) {
        int bound = ((n - 1) * (n - 1)) / 4;
        int y[MAXN], pos[MAXN], arg[MAXN];
        for (int i = 0; i < n; i++) y[i] = i;
        int best = -1; long cnt = 0; long total = 0;
        clock_t t0 = clock();
        do {                                  /* y_0 = 0 fixed: permute the suffix */
            total++;
            for (int i = 0; i < n; i++) pos[y[i]] = i;
            int cur = inv_line(y, n), mn = cur;
            for (int b = 0; b + 1 < n; b++) {
                cur += (n - 1) - 2 * pos[b];
                if (cur < mn) mn = cur;
            }
            if (mn > best) { best = mn; cnt = 1; memcpy(arg, y, sizeof(int) * n); }
            else if (mn == best) cnt++;
        } while (next_perm(y + 1, n - 1));
        double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
        printf("{\"n\":%d,\"bound\":%d,\"line_max\":%d,\"excess\":%d,\"count_y0_0\":%ld,"
               "\"scanned\":%ld,\"arg\":", n, bound, best, best - bound, cnt, total);
        print_perm(arg, n);
        printf(",\"seconds\":%.1f}\n", secs);
        fflush(stdout);
    }
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s perms|lines nmin nmax\n", argv[0]); return 2; }
    int nmin = atoi(argv[2]), nmax = atoi(argv[3]);
    if (nmin < 2 || nmax > MAXN || nmin > nmax) return 2;
    if (!strcmp(argv[1], "perms")) run_perms(nmin, nmax);
    else if (!strcmp(argv[1], "lines")) run_lines(nmin, nmax);
    else return 2;
    return 0;
}
