/* line_profile_fast.c — fast exhaustive computation of I(pi) = min_{q,c}
 * inv(w_{q,c}) for all permutations of n, using O(1)-amortized rotation
 * updates instead of recomputing inv(w) from scratch for each of the n^2
 * cuts (O(n^2) total per permutation instead of the O(n^4) of the reference
 * experiments/line_profile.c).
 *
 * Independent re-implementation for H13-I (session 9), built to extend the
 * exhaustive range of C33 (max_pi I(pi) = floor((n-1)^2/4)) beyond n = 10.
 * Per repo rule 3, cross-checked against experiments/line_profile.c (the
 * trusted reference from session 8) on 4 <= n <= 9: exact match of the full
 * histogram of I(pi) values, including the position of the maximum
 * (checks/check_H13I_session9.py, part 1).
 * Version line_profile_fast-1.0.
 *
 * Algorithm per permutation pi (array of n ints, values 0..n-1):
 *  - For q = 0: v = pi (read starting at position 0, i.e. line for cut q=n-1,
 *    to keep indices simple we instead literally rotate pi itself).
 *  - inv(v) computed once by O(n^2) (n <= 12, fine).
 *  - Sweep c = 0..n-1 using the O(1) update:
 *      j0 = position of value c in v (i.e. v[j0] == c)
 *      new_inv = inv + (n - 1 - 2*j0)
 *      (moving the "wrap point" from value c to value c+1 mod n is exactly
 *      the mirror of rotating positions: the element currently holding the
 *      smallest remaining relabelled value jumps to the largest.)
 *    track running min over c.
 *  - Then rotate v by one step to the left (v <- v[1:] + v[:1], i.e. q += 1),
 *    updating inv(v) by the standard cyclic-rotation formula
 *      inv_new = inv_old + (n - 1 - 2*v[0])
 *    before recomputing the c-sweep for the new q.
 *  - Track the global min over all (q, c) as I(pi).
 *
 * This is O(n) work per q (c-sweep) and O(1) amortized per q-step, so O(n^2)
 * total per permutation instead of O(n^4).
 *
 * Usage: i_fast n [dump.bin]
 * Prints: n, count, max I, one lexicographically-first argmax permutation,
 * and the full histogram of I values.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 13

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
    if (argc < 2) { fprintf(stderr, "usage: %s n [dump.bin]\n", argv[0]); return 2; }
    int n = atoi(argv[1]);
    if (n < 2 || n > MAXN) { fprintf(stderr, "n out of range\n"); return 2; }
    long total = 1;
    for (int k = 2; k <= n; k++) total *= k;

    unsigned char *dump = NULL;
    if (argc >= 3) dump = malloc(total);

    long hist[MAXN * MAXN];
    memset(hist, 0, sizeof hist);
    int global_max = -1;
    long argmax_rank = -1;
    int argmax_perm[MAXN];

    int pi[MAXN], v[MAXN], pos_of_val[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;

    long rank = 0;
    do {
        memcpy(v, pi, n * sizeof(int));
        int inv = 0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                if (v[i] > v[j]) inv++;

        int best = 1 << 30;

        for (int q = 0; q < n; q++) {
            /* c-sweep starting at c = 0: w = v (no shift) has inv = inv(v) */
            for (int i = 0; i < n; i++) pos_of_val[v[i]] = i;
            int cur = inv;
            if (cur < best) best = cur;
            for (int c = 0; c < n - 1; c++) {
                int j0 = pos_of_val[c];
                cur += (n - 1 - 2 * j0);
                if (cur < best) best = cur;
            }
            if (q < n - 1) {
                /* rotate v left by one: q -> q+1 */
                int front = v[0];
                inv += (n - 1 - 2 * front);
                for (int i = 0; i < n - 1; i++) v[i] = v[i + 1];
                v[n - 1] = front;
            }
        }

        if (dump) dump[rank] = (unsigned char)best;
        hist[best]++;
        if (best > global_max) {
            global_max = best;
            argmax_rank = rank;
            memcpy(argmax_perm, pi, n * sizeof(int));
        }
        rank++;
    } while (next_perm(pi, n));

    if (dump) {
        FILE *f = fopen(argv[2], "wb");
        fwrite(dump, 1, total, f);
        fclose(f);
    }

    printf("{\"n\": %d, \"count\": %ld, \"max_I\": %d, \"bound_floor_(n-1)^2/4\": %d, "
           "\"argmax_rank\": %ld, \"argmax_perm\": [", n, total, global_max,
           ((n - 1) * (n - 1)) / 4, argmax_rank);
    for (int i = 0; i < n; i++) printf("%s%d", i ? "," : "", argmax_perm[i]);
    printf("], \"hist\": [");
    int first = 1;
    for (int i = 0; i < n * n; i++) if (hist[i]) {
        printf("%s[%d,%ld]", first ? "" : ",", i, hist[i]);
        first = 0;
    }
    printf("]}\n");
    return 0;
}
