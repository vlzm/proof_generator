/* H13-I, session 9: does 2-dimensional local optimality of the cut grid
 * suffice to bound inv <= floor((n-1)^2/4)?
 *
 * Context: docs/notes/h13_line_model.md lemma C shows that the *necessary*
 * conditions of a GLOBAL minimum (full-axis cumulative sums, one axis at a
 * time) are not sufficient -- there is a linear w at n=7 satisfying them with
 * inv=10 > 9, whose torus class has the true minimum I=8 at a cut that moves
 * BOTH q and c together.  This experiment tests the natural "joint shift"
 * strengthening suggested by that gap: is it enough for (q,c) to be a
 * genuine 2-D local minimum of the discrete grid inv(q,c), i.e. no *single*
 * simultaneous move in q and c (not just axis-aligned) improves it?
 *
 * For every permutation of n, build the full n x n grid inv(q,c) (all double
 * cuts, see experiments/line_model.py for the definition of w and inv), then:
 *
 *  (a) mode "grid": for every cell that is a local minimum under the
 *      4-neighbor (axis-only) and, separately, the 8-neighbor (axis + both
 *      diagonals) torus neighborhoods, record whether inv > bound and track
 *      the worst (max) value found among local minima of each kind.
 *  (b) mode "dist": for every 8-neighbor local minimum that violates the
 *      bound, find the torus L1 distance (in (q,c)) to the nearest cell that
 *      attains the true global minimum I(pi), and track the max and average
 *      such distance -- this measures how "non-local" the true optimum is
 *      relative to a bad local optimum, i.e. whether any FIXED finite radius
 *      could ever repair the argument.
 *
 * Usage: h13_local2d grid n   |   h13_local2d dist n
 * Exhaustive over all n! permutations for n <= 9 (matches C2v table range
 * used elsewhere in the project for exhaustive checks).
 * Version h13_local2d-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 9

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

static int tor(int x, int n) { int a = x < 0 ? -x : x; int b = n - a; return a < b ? a : b; }

static void build_grid(const int *pi, int n, int grid[MAXN][MAXN]) {
    int w[MAXN];
    for (int q = 0; q < n; q++) {
        for (int c = 0; c < n; c++) {
            int shift = q + 1 - c;
            for (int j = 0; j < n; j++) {
                int v = (pi[(q + 1 + j) % n] - shift) % n;
                if (v < 0) v += n;
                w[j] = v;
            }
            grid[q][c] = inversions(w, n);
        }
    }
}

static int run_grid(int n) {
    int bound = ((n - 1) * (n - 1)) / 4;
    static int grid[MAXN][MAXN];
    int pi[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;
    long total = 1;
    for (int k = 2; k <= n; k++) total *= k;

    int worst4 = -1, worst8 = -1;
    long rank = 0, rank_worst4 = -1, rank_worst8 = -1;
    long violations4 = 0, violations8 = 0;
    long total_localmin4 = 0, total_localmin8 = 0;
    int best_pi4[MAXN], best_pi8[MAXN];
    int best_q4 = -1, best_c4 = -1, best_q8 = -1, best_c8 = -1;

    do {
        build_grid(pi, n, grid);
        for (int q = 0; q < n; q++) {
            for (int c = 0; c < n; c++) {
                int v = grid[q][c];
                int qp = (q + 1) % n, qm = (q - 1 + n) % n;
                int cp = (c + 1) % n, cm = (c - 1 + n) % n;
                int is4 = (v <= grid[qp][c] && v <= grid[qm][c] &&
                           v <= grid[q][cp] && v <= grid[q][cm]);
                if (is4) {
                    total_localmin4++;
                    if (v > bound) violations4++;
                    if (v > worst4) {
                        worst4 = v; rank_worst4 = rank; best_q4 = q; best_c4 = c;
                        for (int i = 0; i < n; i++) best_pi4[i] = pi[i];
                    }
                    int is8 = is4 &&
                        v <= grid[qp][cp] && v <= grid[qp][cm] &&
                        v <= grid[qm][cp] && v <= grid[qm][cm];
                    if (is8) {
                        total_localmin8++;
                        if (v > bound) violations8++;
                        if (v > worst8) {
                            worst8 = v; rank_worst8 = rank; best_q8 = q; best_c8 = c;
                            for (int i = 0; i < n; i++) best_pi8[i] = pi[i];
                        }
                    }
                }
            }
        }
        rank++;
    } while (next_perm(pi, n));

    printf("{\"mode\": \"grid\", \"n\": %d, \"bound\": %d, \"perms\": %ld, ", n, bound, total);
    printf("\"local4\": {\"total\": %ld, \"violations\": %ld, \"worst\": %d, \"rank\": %ld, \"q\": %d, \"c\": %d, \"pi\": [",
           total_localmin4, violations4, worst4, rank_worst4, best_q4, best_c4);
    for (int i = 0; i < n; i++) printf("%d%s", best_pi4[i], i + 1 < n ? "," : "");
    printf("]}, ");
    printf("\"local8\": {\"total\": %ld, \"violations\": %ld, \"worst\": %d, \"rank\": %ld, \"q\": %d, \"c\": %d, \"pi\": [",
           total_localmin8, violations8, worst8, rank_worst8, best_q8, best_c8);
    for (int i = 0; i < n; i++) printf("%d%s", best_pi8[i], i + 1 < n ? "," : "");
    printf("]}}\n");
    return 0;
}

static int run_dist(int n) {
    int bound = ((n - 1) * (n - 1)) / 4;
    static int grid[MAXN][MAXN];
    int pi[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;
    long total = 1;
    for (int k = 2; k <= n; k++) total *= k;

    long n_violations = 0;
    int max_dist = -1;
    double sum_dist = 0.0;
    long rank = 0, rank_maxdist = -1;
    int best_pi[MAXN] = {0}; int bq = -1, bc = -1, bgq = -1, bgc = -1;

    do {
        build_grid(pi, n, grid);
        int gmin = 999;
        for (int q = 0; q < n; q++)
            for (int c = 0; c < n; c++)
                if (grid[q][c] < gmin) gmin = grid[q][c];
        for (int q = 0; q < n; q++) {
            for (int c = 0; c < n; c++) {
                int v = grid[q][c];
                if (v <= bound) continue;
                int qp = (q + 1) % n, qm = (q - 1 + n) % n;
                int cp = (c + 1) % n, cm = (c - 1 + n) % n;
                int is8 = (v <= grid[qp][c] && v <= grid[qm][c] &&
                           v <= grid[q][cp] && v <= grid[q][cm] &&
                           v <= grid[qp][cp] && v <= grid[qp][cm] &&
                           v <= grid[qm][cp] && v <= grid[qm][cm]);
                if (!is8) continue;
                n_violations++;
                int best_d = 999, bq2 = -1, bc2 = -1;
                for (int q2 = 0; q2 < n; q2++)
                    for (int c2 = 0; c2 < n; c2++)
                        if (grid[q2][c2] == gmin) {
                            int dd = tor(q2 - q, n) + tor(c2 - c, n);
                            if (dd < best_d) { best_d = dd; bq2 = q2; bc2 = c2; }
                        }
                sum_dist += best_d;
                if (best_d > max_dist) {
                    max_dist = best_d; rank_maxdist = rank;
                    bq = q; bc = c; bgq = bq2; bgc = bc2;
                    for (int i = 0; i < n; i++) best_pi[i] = pi[i];
                }
            }
        }
        rank++;
    } while (next_perm(pi, n));

    printf("{\"mode\": \"dist\", \"n\": %d, \"bound\": %d, \"violations\": %ld, \"max_dist\": %d, "
           "\"avg_dist\": %.3f, \"rank\": %ld, \"q\": %d, \"c\": %d, \"gq\": %d, \"gc\": %d, \"pi\": [",
           n, bound, n_violations, max_dist, n_violations ? sum_dist / n_violations : 0.0,
           rank_maxdist, bq, bc, bgq, bgc);
    for (int i = 0; i < n; i++) printf("%d%s", best_pi[i], i + 1 < n ? "," : "");
    printf("]}\n");
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s grid|dist n\n", argv[0]); return 2; }
    int n = atoi(argv[2]);
    if (n < 2 || n > MAXN) { fprintf(stderr, "n out of range 2..%d\n", MAXN); return 2; }
    if (strcmp(argv[1], "grid") == 0) return run_grid(n);
    if (strcmp(argv[1], "dist") == 0) return run_dist(n);
    fprintf(stderr, "unknown mode %s\n", argv[1]);
    return 2;
}
