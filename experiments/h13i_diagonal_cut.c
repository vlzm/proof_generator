/* h13i_diagonal_cut.c -- exhaustive check of the "fixed c" diagonal family of
 * double cuts for H13-I (docs/notes/h13_line_model.md, session 9).
 *
 * Fix the shift c and let the position cut q range over Z_n (both the position
 * cut and the value cut move together: e = q + 1 - c, a "diagonal" in (q, e)
 * space, i.e. the coupled walk suggested in the session assignment restricted
 * to a straight diagonal instead of an adaptive path):
 *
 *   D_c(pi) = min_{q in Z_n} inv(w_q),   w_q[j] = (pi[(q+1+j) % n] - (q+1-c)) % n
 *
 * Brute force per rotation (O(n^2) per permutation, O(n^2 * n!) total): cheap
 * up to n = 10 (3 628 800 permutations, ~6 s).
 *
 * By the substitution pi -> pi shifted by an additive constant (a bijection of
 * S_n), max_pi D_c(pi) does not depend on c; c = 0 is used below and other c
 * were spot-checked to match (see report).
 *
 * Usage: h13i_diagonal_cut n [c]   (c defaults to 0)
 * Version h13i_diagonal_cut-1.0.
 */
#include <stdio.h>
#include <stdlib.h>

#define MAXN 12

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

static int inversions(const int *w, int n) {
    int c = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (w[i] > w[j]) c++;
    return c;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s n [c]\n", argv[0]); return 2; }
    int n = atoi(argv[1]);
    int c = argc >= 3 ? atoi(argv[2]) : 0;
    if (n < 2 || n > MAXN) return 2;
    int pi[MAXN], w[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;
    long worst = -1;
    int argworst[MAXN];
    long total = 0;
    do {
        total++;
        int best = 999;
        for (int q = 0; q < n; q++) {
            int shift = ((q + 1 - c) % n + n) % n;
            for (int j = 0; j < n; j++) w[j] = ((pi[(q + 1 + j) % n] - shift) % n + n) % n;
            int iv = inversions(w, n);
            if (iv < best) best = iv;
        }
        if (best > worst) { worst = best; for (int i = 0; i < n; i++) argworst[i] = pi[i]; }
    } while (next_perm(pi, n));
    int thresh = (n - 1) * (n - 1) / 4;
    printf("{\"n\": %d, \"c\": %d, \"total\": %ld, \"thresh_floor_(n-1)^2/4\": %d, "
           "\"max_D\": %ld, \"gap\": %ld, \"one_argmax\": [",
           n, c, total, thresh, worst, worst - thresh);
    for (int i = 0; i < n; i++) printf("%d%s", argworst[i], i + 1 < n ? ", " : "");
    printf("]}\n");
    return 0;
}
