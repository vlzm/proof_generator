/* toric_inv.c — H13-I: I(pi) = min over all n^2 double cuts of inv(w),
 * where w_j = pi(a + j) - b (mod n), j = 0..n-1, (a, b) in Z_n x Z_n
 * (a = q + 1, b = q + 1 - c in the (q, c) notation of docs/notes/h13_line_model.md).
 *
 * Two ingredients make an exhaustive scan cheap (both proved in
 * docs/notes/h13i_toric.md, Lemma 1 and Lemma 3):
 *
 *   (1) shift recurrences — all n^2 cuts cost O(n^2) per permutation:
 *         inv(a+1, b) = inv(a, b) + n - 1 - 2 * ((pi(a) - b) mod n)
 *         inv(a, b+1) = inv(a, b) + n - 1 - 2 * ((pi^{-1}(b) - a) mod n)
 *   (2) I is constant on toric classes, so only permutations with pi(0) = 0
 *       are enumerated (exactly one per value rotation: the map
 *       pi -> (pi - pi(0), pi(0)) is a bijection S_n <-> {normalised} x Z_n,
 *       so every count over normalised representatives times n is the count
 *       over all n! permutations).
 *
 * Reported per n: max I, the target floor((n-1)^2/4), the number of
 * normalised permutations attaining the max and the first few of them, the
 * winding number omega(pi) = (1/n) sum_i ((pi(i+1) - pi(i)) mod n) of each
 * maximiser, max I per omega, and the histogram of I.
 *
 * Usage: toric_inv n [--verify] [--json out.json]
 *   --verify  cross-checks the recurrence against a direct O(n^2) inversion
 *             count on every cut of every permutation (use for n <= 8).
 * Version toric_inv-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXN 13
#define MAXI 128

static int n;

static int inversions_direct(const int *w, int m) {
    int inv = 0;
    for (int i = 0; i < m; i++)
        for (int j = i + 1; j < m; j++)
            if (w[i] > w[j]) inv++;
    return inv;
}

/* min over all n^2 cuts, by the shift recurrences */
static int Ival(const int *pi, const int *ipi) {
    int inv00 = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (pi[i] > pi[j]) inv00++;
    int best = inv00;
    int row = inv00;
    for (int a = 0; a < n; a++) {
        if (a > 0) row += n - 1 - 2 * pi[a - 1];
        int cur = row;
        if (cur < best) best = cur;
        for (int b = 0; b + 1 < n; b++) {
            int t = ipi[b] - a;
            if (t < 0) t += n;
            cur += n - 1 - 2 * t;
            if (cur < best) best = cur;
        }
    }
    return best;
}

/* the same by brute force, for --verify */
static int Ival_direct(const int *pi) {
    int w[MAXN];
    int best = 1 << 30;
    for (int a = 0; a < n; a++)
        for (int b = 0; b < n; b++) {
            for (int j = 0; j < n; j++) {
                int v = pi[(a + j) % n] - b;
                v %= n;
                if (v < 0) v += n;
                w[j] = v;
            }
            int inv = inversions_direct(w, n);
            if (inv < best) best = inv;
        }
    return best;
}

static int winding(const int *pi) {
    int s = 0;
    for (int i = 0; i < n; i++) {
        int d = pi[(i + 1) % n] - pi[i];
        d %= n;
        if (d < 0) d += n;
        s += d;
    }
    return s / n;
}

static int next_perm(int *a, int m) {
    int i = m - 2;
    while (i >= 0 && a[i] >= a[i + 1]) i--;
    if (i < 0) return 0;
    int j = m - 1;
    while (a[j] <= a[i]) j--;
    int t = a[i]; a[i] = a[j]; a[j] = t;
    for (int l = i + 1, r = m - 1; l < r; l++, r--) { t = a[l]; a[l] = a[r]; a[r] = t; }
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s n [--verify] [--json out.json]\n", argv[0]); return 2; }
    n = atoi(argv[1]);
    if (n < 3 || n > MAXN) { fprintf(stderr, "n out of range\n"); return 2; }
    int verify = 0;
    const char *jsonpath = NULL;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--verify")) verify = 1;
        else if (!strcmp(argv[i], "--json") && i + 1 < argc) jsonpath = argv[++i];
    }
    int target = ((n - 1) * (n - 1)) / 4;

    int pi[MAXN], ipi[MAXN];
    pi[0] = 0;
    for (int i = 1; i < n; i++) pi[i] = i;

    long total = 0, over = 0, atmax = 0;
    int maxI = -1;
    int maxI_by_omega[MAXN];
    for (int i = 0; i < MAXN; i++) maxI_by_omega[i] = -1;
    long hist[MAXI];
    memset(hist, 0, sizeof hist);
    int argmax[8][MAXN];
    int nargmax = 0;
    int argmax_all_reflection = 1;
    int first_over[MAXN];
    int have_over = 0;

    clock_t t0 = clock();
    do {
        for (int i = 0; i < n; i++) ipi[pi[i]] = i;
        int I = Ival(pi, ipi);
        if (verify) {
            int J = Ival_direct(pi);
            if (I != J) {
                fprintf(stderr, "RECURRENCE MISMATCH at ");
                for (int i = 0; i < n; i++) fprintf(stderr, "%d ", pi[i]);
                fprintf(stderr, ": fast %d direct %d\n", I, J);
                return 1;
            }
        }
        total++;
        if (I < MAXI) hist[I]++;
        int om = winding(pi);
        if (I > maxI_by_omega[om]) maxI_by_omega[om] = I;
        if (I > target) {
            if (!have_over) { memcpy(first_over, pi, sizeof(int) * n); have_over = 1; }
            over++;
        }
        if (I > maxI) {
            maxI = I; atmax = 0; nargmax = 0; argmax_all_reflection = 1;
        }
        if (I == maxI) {
            atmax++;
            if (nargmax < 8) memcpy(argmax[nargmax++], pi, sizeof(int) * n);
            /* reflection test: pi(i) = (h - i) mod n for some h */
            int h = pi[0];
            for (int i = 0; i < n; i++) {
                int e = (h - i) % n; if (e < 0) e += n;
                if (pi[i] != e) { argmax_all_reflection = 0; break; }
            }
        }
        /* next permutation of positions 1..n-1 (pi(0) = 0 stays) */
    } while (next_perm(pi + 1, n - 1));
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;

    printf("n = %d: normalised perms = %ld (= n! / n), max I = %d, target floor((n-1)^2/4) = %d, %s\n",
           n, total, maxI, target, maxI <= target ? "H13-I HOLDS" : "H13-I FAILS");
    printf("  perms with I > target: %ld\n", over);
    if (have_over) {
        printf("  first violator:");
        for (int i = 0; i < n; i++) printf(" %d", first_over[i]);
        printf("\n");
    }
    printf("  normalised perms with I = max: %ld (all perms: %ld); all reflections: %s\n",
           atmax, atmax * n, argmax_all_reflection ? "yes" : "NO");
    for (int k = 0; k < nargmax; k++) {
        printf("  argmax[%d]:", k);
        for (int i = 0; i < n; i++) printf(" %d", argmax[k][i]);
        printf("\n");
    }
    printf("  max I by omega:");
    for (int w = 1; w < n; w++) printf(" %d:%d", w, maxI_by_omega[w]);
    printf("\n  time %.2f s\n", secs);

    if (jsonpath) {
        FILE *f = fopen(jsonpath, "w");
        if (!f) return 1;
        fprintf(f, "{\"version\": \"toric_inv-1.0\", \"n\": %d, \"normalised_perms\": %ld, "
                   "\"all_perms\": %ld, \"max_I\": %d, \"target\": %d, \"holds\": %s, "
                   "\"count_I_gt_target\": %ld, \"count_I_eq_max_normalised\": %ld, "
                   "\"count_I_eq_max_all\": %ld, \"argmax_all_reflections\": %s, \"verify\": %s,\n",
                n, total, total * n, maxI, target, maxI <= target ? "true" : "false",
                over, atmax, atmax * n, argmax_all_reflection ? "true" : "false",
                verify ? "true" : "false");
        fprintf(f, " \"argmax_examples\": [");
        for (int k = 0; k < nargmax; k++) {
            fprintf(f, "%s[", k ? ", " : "");
            for (int i = 0; i < n; i++) fprintf(f, "%s%d", i ? ", " : "", argmax[k][i]);
            fprintf(f, "]");
        }
        fprintf(f, "],\n \"max_I_by_omega\": {");
        for (int w = 1; w < n; w++) fprintf(f, "%s\"%d\": %d", w > 1 ? ", " : "", w, maxI_by_omega[w]);
        fprintf(f, "},\n \"hist_I_normalised\": {");
        int first = 1;
        for (int i = 0; i < MAXI; i++) if (hist[i]) {
            fprintf(f, "%s\"%d\": %ld", first ? "" : ", ", i, hist[i]);
            first = 0;
        }
        fprintf(f, "},\n \"seconds\": %.2f}\n", secs);
        fclose(f);
    }
    return maxI <= target ? 0 : 3;
}
