/* h13i_cut_rules.c — exhaustive tests of candidate proof strategies for
 * H13-I: I(pi) = min_{q,c} inv(w) <= floor((n-1)^2/4) for every pi of Z_n
 * (docs/notes/h13_line_model.md §6). Cut convention identical to
 * experiments/line_profile.c: edge {q, q+1}, value shift = q + 1 - c,
 * w[j] = (pi[(q+1+j) mod n] - shift) mod n.
 *
 * Modes (argv[2]):
 *   ruleA   - "self-referential" cut: for each q, c is fixed so the first
 *             element of the line is relabelled 0 (w_0 = 0); report the
 *             worst case of min_q inv over this restricted n-point family.
 *   rule1d  - drop position freedom (fix q = n-1, natural order), keep only
 *             value-rotation c; report failures against the H13-I threshold
 *             and how many of them are exact reflections pi_h(i) = h-i mod n.
 *   ihist   - full 2D search (both q and c free): histogram of I(pi), used
 *             to confirm the extremal set (exactly n reflections at the
 *             max) and the gap to the next level.
 *
 * Version h13i_cut_rules-1.0. Core: none needed (self-contained inversion
 * count). Usage: h13i_cut_rules n mode
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 12

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

static void make_w(const int *pi, int n, int q, int c, int *w) {
    int shift = q + 1 - c;
    for (int j = 0; j < n; j++) {
        int v = (pi[(q + 1 + j) % n] - shift) % n;
        if (v < 0) v += n;
        w[j] = v;
    }
}

static int is_reflection(const int *pi, int n) {
    for (int h = 0; h < n; h++) {
        int ok = 1;
        for (int i = 0; i < n; i++) {
            int v = ((h - i) % n + n) % n;
            if (pi[i] != v) { ok = 0; break; }
        }
        if (ok) return 1;
    }
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s n {ruleA|rule1d|ihist}\n", argv[0]); return 2; }
    int n = atoi(argv[1]);
    if (n < 2 || n > MAXN) return 2;
    const char *mode = argv[2];
    int threshold = ((n - 1) * (n - 1)) / 4;
    int pi[MAXN], w[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;

    if (strcmp(mode, "ruleA") == 0) {
        long total = 0, fail = 0; int worst = -1; long worst_rank = -1;
        long rank = 0;
        do {
            total++;
            int best = 999999;
            for (int q = 0; q < n; q++) {
                int c = ((q + 1 - pi[(q + 1) % n]) % n + n) % n;
                make_w(pi, n, q, c, w);
                int inv = inversions(w, n);
                if (inv < best) best = inv;
            }
            if (best > threshold) { fail++; if (best > worst) { worst = best; worst_rank = rank; } }
            rank++;
        } while (next_perm(pi, n));
        printf("ruleA n=%d threshold=%d total=%ld fail=%ld worst=%d worst_rank=%ld\n",
               n, threshold, total, fail, worst, worst_rank);
    } else if (strcmp(mode, "rule1d") == 0) {
        long total = 0, fail = 0, fail_refl = 0; int worst = -1; long worst_rank = -1;
        long rank = 0;
        do {
            total++;
            int best = 999999;
            for (int c = 0; c < n; c++) {
                for (int j = 0; j < n; j++) { int v = (pi[j] - c) % n; if (v < 0) v += n; w[j] = v; }
                int inv = inversions(w, n);
                if (inv < best) best = inv;
            }
            if (best > threshold) {
                fail++;
                if (is_reflection(pi, n)) fail_refl++;
                if (best > worst) { worst = best; worst_rank = rank; }
            }
            rank++;
        } while (next_perm(pi, n));
        printf("rule1d n=%d threshold=%d total=%ld fail=%ld fail_reflections=%ld worst=%d worst_rank=%ld\n",
               n, threshold, total, fail, fail_refl, worst, worst_rank);
    } else if (strcmp(mode, "ihist") == 0) {
        long hist[MAXN * MAXN];
        memset(hist, 0, sizeof hist);
        long rank = 0;
        do {
            int best = 999999;
            for (int q = 0; q < n; q++)
                for (int c = 0; c < n; c++) {
                    make_w(pi, n, q, c, w);
                    int inv = inversions(w, n);
                    if (inv < best) best = inv;
                }
            hist[best]++;
            rank++;
        } while (next_perm(pi, n));
        printf("ihist n=%d threshold=%d\n", n, threshold);
        for (int i = 0; i < n * n; i++) if (hist[i]) printf("I=%d count=%ld\n", i, hist[i]);
    } else {
        fprintf(stderr, "unknown mode %s\n", mode);
        return 2;
    }
    return 0;
}
