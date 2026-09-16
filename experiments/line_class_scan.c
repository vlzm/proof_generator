/* line_class_scan.c — exhaustive check of H13-I over toric classes.
 *
 * H13-I: I(pi) = min over all n^2 double cuts (a, b) of inv(a, b) is at most
 * floor((n-1)^2/4) for every permutation pi of Z_n.
 *
 * Two reductions make larger n reachable than in line_profile.c:
 *   (R1) I is constant on the toric class of pi (rotation of positions by
 *        alpha and of values by beta), because inv_{pi'}(a, b) = inv_pi(a-alpha,
 *        b-beta); every class contains a representative with pi(0) = 0, so it is
 *        enough to enumerate the (n-1)! permutations fixing 0.
 *   (R2) the whole n x n table of inv(a, b) is filled with O(1) updates:
 *        inv(a, b+1) = inv(a, b) + (n-1) - 2 * ((pinv[b] - a) mod n),
 *        inv(a+1, b) = inv(a, b) + (n-1) - 2 * ((pi[a] - b) mod n),
 *        so one permutation costs O(n^2) instead of O(n^4).
 *
 * Conventions match experiments/line_profile.c / line_model.py: the cut (a, b)
 * gives the line with position x_i = (i - a) mod n and value y_i = (pi(i) - b)
 * mod n; there a = q+1 and b = q+1-c.
 *
 * Output (stdout, one JSON object): n, number of representatives scanned, the
 * bound floor((n-1)^2/4), max I, how many representatives attain it, how many
 * of those are reflections pi(i) = -i mod n, the lexicographically first
 * maximiser, and the histogram of I.  Exit code 1 if max I exceeds the bound.
 *
 * Mode 1 (third argument) additionally reports the largest inv(a, b) over all
 * cuts that are locally minimal along the four lattice directions of the cut
 * torus, (1,0), (0,1), (1,1), (1,-1) — i.e. inv(a + k*dx, b + k*dy) >= inv(a, b)
 * for every k.  A global minimiser is such a cut, so if that maximum stays at or
 * below the bound, four-direction local minimality already implies H13-I.  The
 * two single directions alone do not (lemma C of docs/notes/h13_line_model.md).
 *
 * Usage: line_class_scan n [checkpoint_every [mode]]
 * Version line_class_scan-1.1.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXN 14
#define MAXI 200

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
    if (argc < 2) { fprintf(stderr, "usage: %s n [checkpoint_every]\n", argv[0]); return 2; }
    int n = atoi(argv[1]);
    if (n < 3 || n > MAXN) { fprintf(stderr, "n out of range\n"); return 2; }
    long long checkpoint = (argc >= 3) ? atoll(argv[2]) : 0;
    int mode = (argc >= 4) ? atoi(argv[3]) : 0;

    int bound = ((n - 1) * (n - 1)) / 4;
    int pi[MAXN], pinv[MAXN];
    static int tab[MAXN][MAXN];
    int lm2 = -1, lm4 = -1;        /* max inv at 2- and 4-direction local minima */
    int lm4_pi[MAXN];
    int lm4_a = -1, lm4_b = -1;
    static const int dirs[4][2] = {{1, 0}, {0, 1}, {1, 1}, {1, -1}};
    long long hist[MAXI];
    memset(hist, 0, sizeof hist);

    /* representatives: pi(0) = 0, pi restricted to 1..n-1 runs over all
     * permutations of {1, ..., n-1} in lexicographic order. */
    pi[0] = 0;
    for (int i = 1; i < n; i++) pi[i] = i;

    int maxI = -1;
    long long count = 0, argcount = 0, refl_count = 0;
    int best[MAXN];
    int first_max_seen = 0;
    long long first_max_rank = -1;
    clock_t t0 = clock();

    for (;;) {
        for (int i = 0; i < n; i++) pinv[pi[i]] = i;

        /* inv(0, 0) directly */
        int inv = 0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                if (pi[i] > pi[j]) inv++;

        int I = inv;
        int cur_a0 = inv;              /* inv(a, 0) as a runs */
        for (int a = 0; a < n; a++) {
            int cur = cur_a0;          /* inv(a, b) as b runs */
            for (int b = 0; b < n; b++) {
                if (cur < I) I = cur;
                if (mode) tab[a][b] = cur;
                int j0 = pinv[b] - a; if (j0 < 0) j0 += n;
                cur += (n - 1) - 2 * j0;
            }
            /* cur is now inv(a, n) = inv(a, 0); step a -> a+1 at b = 0 */
            int w0 = pi[a];            /* (pi(a) - 0) mod n */
            cur_a0 += (n - 1) - 2 * w0;
        }

        if (mode) {
            for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) {
                int v = tab[a][b];
                if (v <= lm2 && v <= lm4) continue;
                int ok[4];
                for (int d = 0; d < 4; d++) {
                    ok[d] = 1;
                    for (int k = 1; k < n && ok[d]; k++) {
                        int aa = (a + k * dirs[d][0]) % n;
                        int bb = ((b + k * dirs[d][1]) % n + n) % n;
                        if (tab[aa][bb] < v) ok[d] = 0;
                    }
                }
                if (ok[0] && ok[1] && v > lm2) lm2 = v;
                if (ok[0] && ok[1] && ok[2] && ok[3] && v > lm4) {
                    lm4 = v; lm4_a = a; lm4_b = b;
                    memcpy(lm4_pi, pi, sizeof(int) * n);
                }
            }
        }

        hist[I]++;
        count++;
        if (I > maxI) { maxI = I; argcount = 0; refl_count = 0; first_max_seen = 0; }
        if (I == maxI) {
            argcount++;
            int is_refl = 1;
            for (int i = 0; i < n; i++) { int r = (n - i) % n; if (pi[i] != r) { is_refl = 0; break; } }
            if (is_refl) refl_count++;
            if (!first_max_seen) {
                memcpy(best, pi, sizeof(int) * n);
                first_max_rank = count - 1;
                first_max_seen = 1;
            }
        }
        if (checkpoint && count % checkpoint == 0) {
            double el = (double)(clock() - t0) / CLOCKS_PER_SEC;
            fprintf(stderr, "checkpoint: %lld reps, max I = %d (bound %d), %.1f s\n",
                    count, maxI, bound, el);
            fflush(stderr);
        }
        if (!next_perm(pi + 1, n - 1)) break;
    }

    double elapsed = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("{\"n\": %d, \"reps\": %lld, \"bound\": %d, \"max_I\": %d, \"holds\": %s,\n",
           n, count, bound, maxI, maxI <= bound ? "true" : "false");
    printf(" \"argmax_reps\": %lld, \"argmax_reflections\": %lld, \"first_argmax_rank\": %lld,\n",
           argcount, refl_count, first_max_rank);
    printf(" \"first_argmax\": [");
    for (int i = 0; i < n; i++) printf("%s%d", i ? ", " : "", best[i]);
    printf("],\n \"seconds\": %.2f,\n", elapsed);
    if (mode) {
        printf(" \"localmin2_max_inv\": %d, \"localmin4_max_inv\": %d,\n", lm2, lm4);
        printf(" \"localmin4_argmax\": {\"a\": %d, \"b\": %d, \"pi\": [", lm4_a, lm4_b);
        for (int i = 0; i < n; i++) printf("%s%d", i ? ", " : "", lm4_pi[i]);
        printf("]},\n");
    }
    printf(" \"hist\": [");
    int first = 1;
    for (int i = 0; i < MAXI; i++) if (hist[i]) {
        printf("%s[%d, %lld]", first ? "" : ", ", i, hist[i]);
        first = 0;
    }
    printf("]}\n");
    return maxI <= bound ? 0 : 1;
}
