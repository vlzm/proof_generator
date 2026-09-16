/* toric_inv.c — exhaustive check of H13-I: I(pi) = min over double cuts (q, c)
 * of inv(w) is at most floor((n-1)^2/4).
 *
 * Speed: the whole n x n grid F(a, b) = inv of the line cut at position a with
 * value origin b is computed in O(n^2) per permutation by the two recurrences
 *   F(a+1, b) = F(a, b) + n - 1 - 2 * ((pi(a) - b) mod n)      (front to back)
 *   F(a, b+1) = F(a, b) + n - 1 - 2 * ((pinv(b) - a) mod n)    (min value to max)
 * instead of O(n^4) by relabelling every line (experiments/line_profile.c).
 * Cut convention of line_profile.c: a = q + 1, b = q + 1 - c, so the grids agree
 * as sets and I(pi) = min over the grid.
 *
 * I(pi) is invariant under the toric action pi -> pi(. + r) - s, so only the
 * (n-1)! permutations with pi(0) = 0 are enumerated (every toric class has at
 * least one such representative).
 *
 * Output: JSON with the histogram of I, the max, the bound, the extremal
 * representatives (up to a cap) and the first violation, if any.
 *
 * Usage: toric_inv n [--full] [--dump-extremal K]
 *   --full            enumerate all n! permutations (cross-check of invariance)
 *   --dump-extremal K print at most K representatives attaining the max
 * Version toric_inv-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 14
#define MAXI 128

static int n;

static int inversions(const int *w, int len) {
    int inv = 0;
    for (int i = 0; i < len; i++)
        for (int j = i + 1; j < len; j++)
            if (w[i] > w[j]) inv++;
    return inv;
}

static int next_perm(int *a, int len) {
    int i = len - 2;
    while (i >= 0 && a[i] >= a[i + 1]) i--;
    if (i < 0) return 0;
    int j = len - 1;
    while (a[j] <= a[i]) j--;
    int t = a[i]; a[i] = a[j]; a[j] = t;
    for (int l = i + 1, r = len - 1; l < r; l++, r--) { t = a[l]; a[l] = a[r]; a[r] = t; }
    return 1;
}

/* minimum of F over the whole grid, and its argmin */
static int grid_min(const int *pi, const int *pinv, int *arga, int *argb) {
    int row[MAXN];
    /* F(0, 0) = inv(pi) */
    int f = inversions(pi, n);
    /* row b = 0..n-1 for a = 0 */
    int best = f, ba = 0, bb = 0;
    row[0] = f;
    for (int b = 1; b < n; b++) {
        f += n - 1 - 2 * pinv[b - 1];          /* (pinv[b-1] - 0) mod n */
        row[b] = f;
        if (f < best) { best = f; ba = 0; bb = b; }
    }
    for (int a = 1; a < n; a++) {
        int v = pi[a - 1];
        for (int b = 0; b < n; b++) {
            int t = v - b; if (t < 0) t += n;
            row[b] += n - 1 - 2 * t;
            if (row[b] < best) { best = row[b]; ba = a; bb = b; }
        }
    }
    *arga = ba; *argb = bb;
    return best;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s n [--full] [--dump-extremal K]\n", argv[0]); return 2; }
    n = atoi(argv[1]);
    if (n < 3 || n > MAXN) { fprintf(stderr, "n out of range\n"); return 2; }
    int full = 0, dumpk = 5;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--full")) full = 1;
        else if (!strcmp(argv[i], "--dump-extremal") && i + 1 < argc) dumpk = atoi(argv[++i]);
    }
    int bound = ((n - 1) * (n - 1)) / 4;

    int pi[MAXN], pinv[MAXN];
    long hist[MAXI]; memset(hist, 0, sizeof hist);
    int maxI = -1; long nmax = 0, total = 0, viol = 0;
    int ex[64][MAXN]; int nex = 0;
    int firstviol[MAXN]; int haveviol = 0;

    /* enumerate: full -> all n!, else pi(0) = 0 and the rest permuted */
    int start = full ? 0 : 1;
    for (int i = 0; i < n; i++) pi[i] = i;
    do {
        for (int i = 0; i < n; i++) pinv[pi[i]] = i;
        int aa, bb;
        int I = grid_min(pi, pinv, &aa, &bb);
        total++;
        if (I >= 0 && I < MAXI) hist[I]++;
        if (I > maxI) { maxI = I; nmax = 0; nex = 0; }
        if (I == maxI) {
            nmax++;
            if (nex < dumpk && nex < 64) { memcpy(ex[nex], pi, sizeof(int) * n); nex++; }
        }
        if (I > bound) {
            viol++;
            if (!haveviol) { memcpy(firstviol, pi, sizeof(int) * n); haveviol = 1; }
        }
    } while (next_perm(pi + start, n - start));

    printf("{\"version\": \"toric_inv-1.0\", \"n\": %d, \"enumerated\": %ld, \"scope\": \"%s\",\n",
           n, total, full ? "all n! permutations" : "pi(0)=0 representatives of toric classes");
    printf(" \"bound_floor_(n-1)^2/4\": %d, \"max_I\": %d, \"count_max\": %ld, \"violations\": %ld,\n",
           bound, maxI, nmax, viol);
    if (haveviol) {
        printf(" \"first_violation\": [");
        for (int i = 0; i < n; i++) printf("%s%d", i ? ", " : "", firstviol[i]);
        printf("],\n");
    }
    printf(" \"extremal_examples\": [");
    for (int k = 0; k < nex; k++) {
        printf("%s[", k ? ", " : "");
        for (int i = 0; i < n; i++) printf("%s%d", i ? ", " : "", ex[k][i]);
        printf("]");
    }
    printf("],\n \"hist_I\": [");
    int first = 1;
    for (int i = 0; i < MAXI; i++) if (hist[i]) { printf("%s[%d, %ld]", first ? "" : ", ", i, hist[i]); first = 0; }
    printf("]}\n");
    return 0;
}
