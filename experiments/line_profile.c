/* line_profile.c — I(pi) = min over cut q and shift c of the inversion count of
 * the relabelled line (see experiments/line_model.py), for every permutation
 * of n (Lehmer / lexicographic rank order), joined with the distance table
 * data/tables/dist_n{n}.bin.
 *
 * Output (stdout, one JSON object): for each value I: count of pi, max d, min d,
 * the lexicographically first argmax; global max of d - 2I and d - I with
 * argmax; and the reverse profile: for each d, max and min I.
 *
 * Usage: line_profile n path/to/dist_n{n}.bin [path/to/I_n{n}.bin]
 * The optional third argument dumps I(pi) as one byte per permutation.
 * Version line_profile-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 12
#define MAXI 80
#define MAXD 80

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

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s n dist.bin [I.bin]\n", argv[0]); return 2; }
    int n = atoi(argv[1]);
    if (n < 2 || n > MAXN) return 2;
    long total = 1;
    for (int k = 2; k <= n; k++) total *= k;
    unsigned char *dist = malloc(total);
    FILE *f = fopen(argv[2], "rb");
    if (!f || fread(dist, 1, total, f) != (size_t)total) { fprintf(stderr, "cannot read %s\n", argv[2]); return 1; }
    fclose(f);
    unsigned char *Ivals = NULL;
    if (argc >= 4) Ivals = malloc(total);

    long cnt[MAXI]; int maxd[MAXI], mind[MAXI]; long argmax[MAXI];
    int maxI_of_d[MAXD], minI_of_d[MAXD]; long cnt_d[MAXD];
    for (int i = 0; i < MAXI; i++) { cnt[i] = 0; maxd[i] = -1; mind[i] = 999; argmax[i] = -1; }
    for (int i = 0; i < MAXD; i++) { maxI_of_d[i] = -1; minI_of_d[i] = 999; cnt_d[i] = 0; }
    int best_ex2 = -999; long arg_ex2 = -1; int best_ex1 = -999; long arg_ex1 = -1;
    /* joint histogram of (I, d) */
    static long joint[MAXI][MAXD];
    memset(joint, 0, sizeof joint);

    int pi[MAXN], w[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;
    long rank = 0;
    do {
        int best = 999;
        for (int q = 0; q < n; q++) {
            for (int c = 0; c < n; c++) {
                int shift = q + 1 - c;
                for (int j = 0; j < n; j++) {
                    int v = (pi[(q + 1 + j) % n] - shift) % n;
                    if (v < 0) v += n;
                    w[j] = v;
                }
                int inv = inversions(w, n);
                if (inv < best) best = inv;
            }
        }
        int d = dist[rank];
        if (Ivals) Ivals[rank] = (unsigned char)best;
        cnt[best]++;
        joint[best][d]++;
        cnt_d[d]++;
        if (d > maxd[best]) { maxd[best] = d; argmax[best] = rank; }
        if (d < mind[best]) mind[best] = d;
        if (best > maxI_of_d[d]) maxI_of_d[d] = best;
        if (best < minI_of_d[d]) minI_of_d[d] = best;
        if (d - 2 * best > best_ex2) { best_ex2 = d - 2 * best; arg_ex2 = rank; }
        if (d - best > best_ex1) { best_ex1 = d - best; arg_ex1 = rank; }
        rank++;
    } while (next_perm(pi, n));

    if (Ivals) { f = fopen(argv[3], "wb"); fwrite(Ivals, 1, total, f); fclose(f); }

    printf("{\"n\": %d, \"count\": %ld, \"max_d_minus_2I\": %d, \"argmax_d_minus_2I_rank\": %ld, "
           "\"max_d_minus_I\": %d, \"argmax_d_minus_I_rank\": %ld,\n", n, total, best_ex2, arg_ex2, best_ex1, arg_ex1);
    printf(" \"by_I\": [");
    int first = 1;
    for (int i = 0; i < MAXI; i++) if (cnt[i]) {
        printf("%s\n  {\"I\": %d, \"count\": %ld, \"max_d\": %d, \"min_d\": %d, \"argmax_rank\": %ld}", first ? "" : ",", i, cnt[i], maxd[i], mind[i], argmax[i]);
        first = 0;
    }
    printf("],\n \"by_d\": [");
    first = 1;
    for (int d = 0; d < MAXD; d++) if (cnt_d[d]) {
        printf("%s\n  {\"d\": %d, \"count\": %ld, \"max_I\": %d, \"min_I\": %d}", first ? "" : ",", d, cnt_d[d], maxI_of_d[d], minI_of_d[d]);
        first = 0;
    }
    printf("],\n \"joint\": [");
    first = 1;
    for (int i = 0; i < MAXI; i++) for (int d = 0; d < MAXD; d++) if (joint[i][d]) {
        printf("%s[%d, %d, %ld]", first ? "" : ", ", i, d, joint[i][d]);
        first = 0;
    }
    printf("]}\n");
    free(dist);
    return 0;
}
