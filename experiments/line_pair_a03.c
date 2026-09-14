/* line_pair_a03.c -- H13-I, session 9: exhaustive check of the fixed
 * position-cut pair {a1, a2} (default 0, 3), independent of pi: does
 * min(min_b inv(a1,b), min_b inv(a2,b)) <= floor((n-1)^2/4) for every
 * permutation of n? Companion to experiments/h13_averaging_bound.py (which
 * covers this exhaustively only up to n = 8 in pure Python; this program
 * extends the exhaustive range to n = 9..11, where a Python loop is too slow).
 * The pair {0, 3} is NOT claimed to work in general -- see the hill-climbing
 * counterexample at n = 20 in h13_averaging_bound.py (excess 2): this program
 * only certifies (or refutes) it exhaustively for the given n.
 * Usage: line_pair_a03 n a1 a2
 * Version line_pair_a03-1.0.
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
static int minb(const int *pi, int n, int a) {
    int best = 99999, w[MAXN];
    for (int b = 0; b < n; b++) {
        for (int j = 0; j < n; j++) { int v = (pi[(a + j) % n] - b) % n; if (v < 0) v += n; w[j] = v; }
        int iv = inversions(w, n);
        if (iv < best) best = iv;
    }
    return best;
}
int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s n a1 a2\n", argv[0]); return 2; }
    int n = atoi(argv[1]), a1 = atoi(argv[2]), a2 = atoi(argv[3]);
    if (n < 2 || n > MAXN) return 2;
    int pi[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;
    int thr = (n - 1) * (n - 1) / 4;
    int maxV = -1; long cnt = 0, fails = 0, total = 0;
    int argmax[MAXN];
    do {
        int v1 = minb(pi, n, a1);
        int v2 = (a2 == a1) ? v1 : minb(pi, n, a2);
        int m = v1 < v2 ? v1 : v2;
        if (m > thr) fails++;
        if (m > maxV) { maxV = m; cnt = 1; for (int i = 0; i < n; i++) argmax[i] = pi[i]; }
        else if (m == maxV) cnt++;
        total++;
    } while (next_perm(pi, n));
    printf("{\"n\": %d, \"a1\": %d, \"a2\": %d, \"threshold\": %d, \"worst\": %d, \"excess\": %d, "
           "\"fails\": %ld, \"cnt_at_worst\": %ld, \"total\": %ld, \"argworst\": [",
           n, a1, a2, thr, maxV, maxV - thr, fails, cnt, total);
    for (int i = 0; i < n; i++) printf("%d%s", argmax[i], i + 1 < n ? "," : "");
    printf("]}\n");
    return 0;
}
