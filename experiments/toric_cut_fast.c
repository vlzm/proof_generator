/* toric_cut_fast.c — O(n^2)-per-permutation computation of
 * I(pi) = min_{q,c} inv(w_{q,c}), the double-cut inversion minimum used in
 * H13-I (docs/notes/h13_line_model.md §6, PLAN.md H13). Replaces the
 * O(n^4)-per-permutation loop of experiments/line_profile.c (all n^2 cuts,
 * each inversion count O(n^2)) with an O(n^2) incremental recursion, derived
 * and proved in docs/proofs/D1_toric_cut_recursion.md.
 *
 * Convention (equivalent to line_model.py's (q,c) up to relabelling,
 * checked against the brute-force line_profile algorithm below for
 * 4 <= n <= 8, exhaustively, mode "check"):
 *   w^{(a,b)}_i = (pi[(a+i) mod n] - b) mod n,  i = 0..n-1
 *   Inv(a,b) = inversions of w^{(a,b)} as a linear sequence
 *   I(pi) = min_{a,b in Z_n} Inv(a,b)
 *
 * Recursions (proved in D1):
 *   Inv(a+1,b) = Inv(a,b) + (n-1-2*((pi[a]-b) mod n))
 *   Inv(a,b+1) = Inv(a,b) + (n-1-2*((inv_pi[b]-a) mod n))
 *
 * Modes:
 *   check n        — brute force vs fast, all permutations, abort on mismatch
 *   profile n      — exhaustive I(pi) histogram + max I argmax
 * Usage: toric_cut_fast check|profile n
 * Version toric_cut_fast-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 14

static int brute_best(const int *pi, int n) {
    int best = 1000000;
    for (int a = 0; a < n; a++) {
        for (int b = 0; b < n; b++) {
            int w[MAXN];
            for (int i = 0; i < n; i++) {
                int v = (pi[(a + i) % n] - b) % n;
                if (v < 0) v += n;
                w[i] = v;
            }
            int inv = 0;
            for (int i = 0; i < n; i++)
                for (int j = i + 1; j < n; j++)
                    if (w[i] > w[j]) inv++;
            if (inv < best) best = inv;
        }
    }
    return best;
}

static int fast_best(const int *pi, const int *inv_pi, int n) {
    /* Inv(0,0) by brute O(n^2) once. */
    int inv00 = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (pi[i] > pi[j]) inv00++;

    int Aa[MAXN]; /* Inv(a,0) */
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
    if (argc < 3) { fprintf(stderr, "usage: %s check|profile n\n", argv[0]); return 2; }
    const char *mode = argv[1];
    int n = atoi(argv[2]);
    if (n < 2 || n > MAXN) { fprintf(stderr, "n out of range\n"); return 2; }

    int pi[MAXN], inv_pi[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;

    if (strcmp(mode, "check") == 0) {
        long total = 0, mismatches = 0;
        do {
            for (int i = 0; i < n; i++) inv_pi[pi[i]] = i;
            int b1 = brute_best(pi, n);
            int b2 = fast_best(pi, inv_pi, n);
            total++;
            if (b1 != b2) {
                mismatches++;
                if (mismatches <= 5) {
                    printf("MISMATCH pi=(");
                    for (int i = 0; i < n; i++) printf("%d%s", pi[i], i + 1 < n ? "," : "");
                    printf(") brute=%d fast=%d\n", b1, b2);
                }
            }
        } while (next_perm(pi, n));
        printf("{\"n\": %d, \"total\": %ld, \"mismatches\": %ld}\n", n, total, mismatches);
        return mismatches ? 1 : 0;
    } else if (strcmp(mode, "profile") == 0) {
        long total = 0;
        int maxI = -1;
        long argmax = -1;
        long cnt[MAXN * MAXN];
        memset(cnt, 0, sizeof cnt);
        long rank = 0;
        do {
            for (int i = 0; i < n; i++) inv_pi[pi[i]] = i;
            int I = fast_best(pi, inv_pi, n);
            cnt[I]++;
            if (I > maxI) { maxI = I; argmax = rank; }
            total++;
            rank++;
        } while (next_perm(pi, n));
        printf("{\"n\": %d, \"total\": %ld, \"maxI\": %d, \"floor_n1_sq_4\": %d, \"argmax_rank\": %ld, \"by_I\": [",
               n, total, maxI, ((n - 1) * (n - 1)) / 4, argmax);
        for (int i = 0; i <= maxI; i++) if (cnt[i]) printf("[%d,%ld],", i, cnt[i]);
        printf("null]}\n");
        return 0;
    }
    fprintf(stderr, "unknown mode\n");
    return 2;
}
