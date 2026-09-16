/* toric_degree.c — inversion profile of the toric class of a permutation and
 * the inversion degrees of its elements at the optimal cuts.
 *
 * For a permutation pi of Z_n and a cut (a, b) the line is
 *     w_j = (pi[(a + j) mod n] - b) mod n,  j = 0..n-1,
 * I(pi) = min over the n^2 cuts of inv(w)  (see docs/notes/h13_line_model.md).
 * For a line w the inversion degree of position j is
 *     deg_j = #{k : the pair (j, k) is inverted},
 * and h_j = (n - 1) - 2 deg_j, so inv(w) = (1/4) (n(n-1) - sum_j h_j).
 *
 * Checked here, exhaustively over all permutations (or over a random sample):
 *   H13-I   I(pi) <= floor((n-1)^2/4);
 *   L       at every cut attaining I(pi): max_j deg_j <= floor((n-1)/2)
 *           (equivalently h_j >= 0 for all j).
 * L implies H13-I for even n, since then h_j is odd, hence h_j >= 1.
 *
 * Also measured (averaging obstruction):
 *   Sigma(pi) = sum over the n^2 cuts of inv  =  n^2 C(n,2) - S(pi),
 *   S(pi) = sum over ordered pairs x != y of ((y-x) mod n)*((pi y - pi x) mod n);
 *   plain averaging proves H13-I for pi iff Sigma(pi) <= n^2 floor((n-1)^2/4).
 *   Claim C37: S(pi) >= n^2 (n^2-1)/6 with equality exactly on the reflections.
 *
 * Usage: toric_degree n [sample_count seed]
 * Output: one JSON object on stdout.  Version toric_degree-1.1.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXN 13

static int n;
static int perm[MAXN], pos[MAXN], line[MAXN], deg[MAXN];
static int profile[MAXN * MAXN];

static int inversions(const int *w) {
    int c = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (w[i] > w[j]) c++;
    return c;
}

/* fill profile[a*n+b] = inv of the line at cut (a, b) */
static void fill_profile(void) {
    int w[MAXN];
    for (int i = 0; i < n; i++) w[i] = perm[i];
    int inv0 = inversions(w);
    for (int a = 0; a < n; a++) {
        if (a > 0) inv0 += (n - 1) - 2 * perm[a - 1];
        int cur = inv0;
        profile[a * n + 0] = cur;
        for (int b = 1; b < n; b++) {
            int j = pos[b - 1] - a;
            if (j < 0) j += n;
            cur += (n - 1) - 2 * j;
            profile[a * n + b] = cur;
        }
    }
}

static void build_line(int a, int b) {
    for (int j = 0; j < n; j++) {
        int v = perm[(a + j) % n] - b;
        if (v < 0) v += n;
        line[j] = v;
    }
}

static void degrees(void) {
    memset(deg, 0, sizeof(int) * n);
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (line[i] > line[j]) { deg[i]++; deg[j]++; }
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

static unsigned long long rng_state;
static unsigned long long rng_next(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s n [sample seed]\n", argv[0]); return 2; }
    n = atoi(argv[1]);
    if (n < 3 || n > MAXN) { fprintf(stderr, "n out of range\n"); return 2; }
    long long sample = 0;
    if (argc >= 3) sample = atoll(argv[2]);
    rng_state = (argc >= 4) ? strtoull(argv[3], NULL, 10) : 20260916ULL;
    if (rng_state == 0) rng_state = 1;

    const int bound = (n - 1) * (n - 1) / 4;      /* H13-I bound */
    const int degbound = (n - 1) / 2;             /* L bound */

    int maxI = -1, worst_deg_all = -1, worst_deg_best = -1;
    long long nperm = 0, nmaxI = 0;
    long long minS = -1, nminS = 0, avg_ok = 0;
    long long S_refl = (long long)n * n * ((long long)n * n - 1) / 6;
    int ex_minS[MAXN]; int have_minS = 0;
    long long viol_H = 0, viol_L_all = 0, viol_L_some = 0;
    int ex_H[MAXN], ex_La[MAXN], ex_Ls[MAXN];
    int have_H = 0, have_La = 0, have_Ls = 0;
    int argmaxI[MAXN]; int have_argmaxI = 0;

    clock_t t0 = clock();
    for (int i = 0; i < n; i++) perm[i] = i;
    for (;;) {
        if (sample > 0) {
            /* random permutation (Fisher-Yates) */
            for (int i = 0; i < n; i++) perm[i] = i;
            for (int i = n - 1; i > 0; i--) {
                int j = (int)(rng_next() % (unsigned long long)(i + 1));
                int t = perm[i]; perm[i] = perm[j]; perm[j] = t;
            }
        }
        for (int i = 0; i < n; i++) pos[perm[i]] = i;
        fill_profile();
        int I = profile[0];
        long long sigma = 0;
        for (int k = 0; k < n * n; k++) { if (profile[k] < I) I = profile[k]; sigma += profile[k]; }
        long long Sval = (long long)n * n * ((long long)n * (n - 1) / 2) - sigma;
        if (minS < 0 || Sval < minS) { minS = Sval; nminS = 0; memcpy(ex_minS, perm, sizeof(int) * n); have_minS = 1; }
        if (Sval == minS) nminS++;
        if (sigma <= (long long)n * n * bound) avg_ok++;
        nperm++;
        if (I > maxI) { maxI = I; nmaxI = 0; memcpy(argmaxI, perm, sizeof(int) * n); have_argmaxI = 1; }
        if (I == maxI) nmaxI++;
        if (I > bound) {
            viol_H++;
            if (!have_H) { memcpy(ex_H, perm, sizeof(int) * n); have_H = 1; }
        }
        /* inversion degrees at the optimal cuts */
        int bad_all = 0, ok_some = 0, best_here = 1 << 30;
        for (int a = 0; a < n; a++)
            for (int b = 0; b < n; b++) {
                if (profile[a * n + b] != I) continue;
                build_line(a, b);
                degrees();
                int md = 0;
                for (int j = 0; j < n; j++) if (deg[j] > md) md = deg[j];
                if (md > worst_deg_all) worst_deg_all = md;
                if (md < best_here) best_here = md;
                if (md > degbound) bad_all = 1; else ok_some = 1;
            }
        if (best_here > worst_deg_best) worst_deg_best = best_here;
        if (bad_all) {
            viol_L_all++;
            if (!have_La) { memcpy(ex_La, perm, sizeof(int) * n); have_La = 1; }
        }
        if (!ok_some) {
            viol_L_some++;
            if (!have_Ls) { memcpy(ex_Ls, perm, sizeof(int) * n); have_Ls = 1; }
        }
        if (sample > 0) { if (nperm >= sample) break; }
        else if (!next_perm(perm)) break;
    }
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;

    printf("{\"version\": \"toric_degree-1.1\", \"n\": %d, \"mode\": \"%s\", \"seed\": %s,\n",
           n, sample > 0 ? "sample" : "exhaustive", argc >= 4 ? argv[3] : "20260916");
    printf(" \"permutations\": %lld, \"bound_H13I\": %d, \"max_I\": %d, \"count_max_I\": %lld,\n",
           nperm, bound, maxI, nmaxI);
    printf(" \"violations_H13I\": %lld, \"deg_bound_L\": %d,\n", viol_H, degbound);
    printf(" \"max_deg_over_all_optimal_cuts\": %d, \"max_over_pi_of_min_over_optimal_cuts_deg\": %d,\n",
           worst_deg_all, worst_deg_best);
    printf(" \"violations_L_at_some_optimal_cut\": %lld, \"violations_L_at_every_optimal_cut\": %lld,\n",
           viol_L_all, viol_L_some);
    printf(" \"min_S\": %lld, \"S_reflection_n2_n2m1_over6\": %lld, \"count_min_S\": %lld,\n",
           minS, S_refl, nminS);
    printf(" \"averaging_suffices_count\": %lld,\n", avg_ok);
    if (have_minS) {
        printf(" \"argmin_S\": [");
        for (int i = 0; i < n; i++) printf("%s%d", i ? ", " : "", ex_minS[i]);
        printf("],\n");
    }
    if (have_argmaxI) {
        printf(" \"argmax_I\": [");
        for (int i = 0; i < n; i++) printf("%s%d", i ? ", " : "", argmaxI[i]);
        printf("],\n");
    }
    if (have_H) {
        printf(" \"example_violating_H13I\": [");
        for (int i = 0; i < n; i++) printf("%s%d", i ? ", " : "", ex_H[i]);
        printf("],\n");
    }
    if (have_La) {
        printf(" \"example_L_fails_at_some_optimal_cut\": [");
        for (int i = 0; i < n; i++) printf("%s%d", i ? ", " : "", ex_La[i]);
        printf("],\n");
    }
    if (have_Ls) {
        printf(" \"example_L_fails_at_every_optimal_cut\": [");
        for (int i = 0; i < n; i++) printf("%s%d", i ? ", " : "", ex_Ls[i]);
        printf("],\n");
    }
    printf(" \"seconds\": %.2f}\n", secs);
    return 0;
}
