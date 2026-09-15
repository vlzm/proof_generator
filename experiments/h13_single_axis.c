/* h13_single_axis.c — checks whether a single-parameter (one circular
 * rotation, n choices) reduction of the double cut (q, c), n^2 choices,
 * see experiments/line_model.py / line_profile.c) already attains
 * H13-I's target bound floor((n-1)^2/4) on max_pi min inv.
 *
 * Two single-axis families, both special cases of the general double cut
 * w_j = pi(q+1+j) - (q+1-c) mod n:
 *   rotate: c = q+1 (no value shift/mod correction) -- w is just a cyclic
 *     rotation of pi read as a plain sequence (no wraparound subtraction).
 *   shift:  q = n-1 (no position rotation) -- w_j = (pi(j) + c) mod n,
 *     c ranges over all n values (value axis rotated, wraparound applied).
 * Exhaustive over all permutations of {0,...,n-1} (lexicographic order).
 *
 * Usage: h13_single_axis n [rotate|shift|both]
 * Version h13_single_axis-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 13

static int n;
static int perm[MAXN];
static int fen[MAXN + 1];

static void fen_reset(int m) { for (int i = 0; i <= m; i++) fen[i] = 0; }
static void fen_add(int i, int v, int m) { for (i++; i <= m; i += i & (-i)) fen[i] += v; }
static int fen_sum(int i) { int s = 0; for (i++; i > 0; i -= i & (-i)) s += fen[i]; return s; }

static int inversions_direct(const int *a) {
    int c = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (a[i] > a[j]) c++;
    return c;
}

static long long best_rotate = -1;
static int worst_rotate[MAXN];
static long long best_shift = -1;
static int worst_shift[MAXN];

static void process_rotate(void) {
    fen_reset(n);
    long long inv = 0;
    for (int i = n - 1; i >= 0; i--) { inv += fen_sum(perm[i] - 1); fen_add(perm[i], 1, n); }
    long long minv = inv;
    for (int k = 0; k < n - 1; k++) {
        inv = inv + (n - 1) - 2LL * perm[k];
        if (inv < minv) minv = inv;
    }
    if (minv > best_rotate) { best_rotate = minv; memcpy(worst_rotate, perm, sizeof(int) * n); }
}

static void process_shift(void) {
    int w[MAXN];
    long long minv = -1;
    for (int c = 0; c < n; c++) {
        for (int j = 0; j < n; j++) w[j] = (perm[j] + c) % n;
        int iv = inversions_direct(w);
        if (minv < 0 || iv < minv) minv = iv;
    }
    if (minv > best_shift) { best_shift = minv; memcpy(worst_shift, perm, sizeof(int) * n); }
}

static int next_perm(int *a) {
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
    if (argc < 2) { fprintf(stderr, "usage: %s n [rotate|shift|both]\n", argv[0]); return 2; }
    n = atoi(argv[1]);
    if (n < 2 || n > MAXN) { fprintf(stderr, "n out of range\n"); return 2; }
    const char *mode = argc >= 3 ? argv[2] : "both";
    int do_rotate = strcmp(mode, "shift") != 0;
    int do_shift = strcmp(mode, "rotate") != 0;

    for (int i = 0; i < n; i++) perm[i] = i;
    do {
        if (do_rotate) process_rotate();
        if (do_shift) process_shift();
    } while (next_perm(perm));

    long long target = (long long)(n - 1) * (n - 1) / 4;
    if (do_rotate) {
        printf("rotate: n=%d max_pi min_k inv = %lld  target floor((n-1)^2/4)=%lld  %s\n",
               n, best_rotate, target, best_rotate <= target ? "OK" : "VIOLATION");
        printf("  worst: "); for (int i = 0; i < n; i++) printf("%d ", worst_rotate[i]); printf("\n");
    }
    if (do_shift) {
        printf("shift:  n=%d max_pi min_c inv = %lld  target floor((n-1)^2/4)=%lld  %s\n",
               n, best_shift, target, best_shift <= target ? "OK" : "VIOLATION");
        printf("  worst: "); for (int i = 0; i < n; i++) printf("%d ", worst_shift[i]); printf("\n");
    }
    return 0;
}
