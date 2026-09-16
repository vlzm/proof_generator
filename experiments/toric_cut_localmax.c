/* toric_cut_localmax.c — tests whether I(pi) has any strict local maximum
 * (under adjacent value-swaps and/or adjacent position-swaps) below the
 * global maximum floor((n-1)^2/4), for H13-I (docs/notes/h13_line_model.md
 * §6). A "no non-global local max" result on all pi would support a
 * compression/exchange proof of H13-I (reflections are the unique fixed
 * points of an I-non-decreasing move, hence the global maximum).
 *
 * For every permutation pi (exhaustive), and for every one of its (n-1)
 * adjacent-VALUE-swap neighbours (swap values v, v+1 for v=0..n-2) and,
 * if requested, (n-1) adjacent-POSITION-swap neighbours (swap positions
 * i, i+1 for i=0..n-2), compute I(neighbour) using the O(n^2) recursion
 * (D1) and report: is there any pi with I(pi) < floor((n-1)^2/4) all of
 * whose checked neighbours have I(neighbour) <= I(pi) (a non-global local
 * max)?
 *
 * Usage: toric_cut_localmax n [value|position|both]
 * Version toric_cut_localmax-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 14

static int inv_of(const int *w, int n) {
    int inv = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (w[i] > w[j]) inv++;
    return inv;
}

static int fast_best(const int *pi, int n) {
    int inv_pi[MAXN];
    for (int i = 0; i < n; i++) inv_pi[pi[i]] = i;
    int inv00 = inv_of(pi, n);
    int Aa[MAXN];
    Aa[0] = inv00;
    for (int a = 0; a < n - 1; a++) {
        int val = pi[a] % n;
        Aa[a + 1] = Aa[a] + (n - 1 - 2 * val);
    }
    int best = Aa[0];
    for (int a = 0; a < n; a++) {
        if (Aa[a] < best) best = Aa[a];
        int cur = Aa[a];
        for (int b = 0; b < n - 1; b++) {
            int istar = (inv_pi[b] - a) % n;
            if (istar < 0) istar += n;
            cur = cur + (n - 1 - 2 * istar);
            if (cur < best) best = cur;
        }
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
    if (argc < 3) { fprintf(stderr, "usage: %s n value|position|both\n", argv[0]); return 2; }
    int n = atoi(argv[1]);
    if (n < 2 || n > MAXN) return 2;
    int mode_val = strcmp(argv[2], "value") == 0 || strcmp(argv[2], "both") == 0;
    int mode_pos = strcmp(argv[2], "position") == 0 || strcmp(argv[2], "both") == 0;
    int bound = ((n - 1) * (n - 1)) / 4;

    int pi[MAXN], q[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;

    long total = 0, local_max_nonglobal = 0;
    do {
        int I0 = fast_best(pi, n);
        total++;
        if (I0 >= bound) continue; /* only interested in sub-maximal pi */
        int improved = 0;
        if (mode_val) {
            for (int v = 0; v < n - 1 && !improved; v++) {
                memcpy(q, pi, sizeof(int) * n);
                for (int i = 0; i < n; i++) {
                    if (q[i] == v) q[i] = v + 1;
                    else if (q[i] == v + 1) q[i] = v;
                }
                if (fast_best(q, n) > I0) improved = 1;
            }
        }
        if (mode_pos && !improved) {
            for (int i = 0; i < n - 1 && !improved; i++) {
                memcpy(q, pi, sizeof(int) * n);
                int t = q[i]; q[i] = q[i + 1]; q[i + 1] = t;
                if (fast_best(q, n) > I0) improved = 1;
            }
        }
        if (!improved) {
            local_max_nonglobal++;
            if (local_max_nonglobal <= 8) {
                printf("NONGLOBAL_LOCAL_MAX pi=(");
                for (int i = 0; i < n; i++) printf("%d%s", pi[i], i + 1 < n ? "," : "");
                printf(") I=%d bound=%d\n", I0, bound);
            }
        }
    } while (next_perm(pi, n));
    printf("{\"n\": %d, \"mode\": \"%s\", \"total\": %ld, \"nonglobal_local_max\": %ld, \"bound\": %d}\n",
           n, argv[2], total, local_max_nonglobal, bound);
    return 0;
}
