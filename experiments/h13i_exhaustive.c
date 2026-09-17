/* h13i_exhaustive.c — exhaustive check of H13-I: I(pi) <= floor((n-1)^2/4)
 * for every permutation pi of {0,...,n-1}, where
 *   I(pi) = min_{a,b in Z_n} inv( (pi(a+j) - b) mod n , j = 0..n-1 )
 * (double toric cut: rotate positions by a, rotate values by b, count
 * ordinary linear inversions of the resulting sequence).
 *
 * Uses the double-incremental recurrence (verified against brute force in
 * scratch/h13i_fast2.py, matches on n = 3..8):
 *   - rotating positions by 1 (move front element to back) changes inv by
 *     (n - 1 - 2*w_0), where w_0 is the value currently at the front;
 *   - rotating values by 1 (subtract 1 mod n from every value) changes inv
 *     by (n - 1 - 2*j*), where j* is the position of the value that is
 *     about to become 0 (i.e. currently equal to the old shift amount b).
 * This gives O(n) work per (a,b) sweep after one O(n^2) baseline inversion
 * count per permutation, i.e. O(n^2) total per permutation instead of the
 * O(n^4) naive re-count over all n^2 cuts (see experiments/line_profile.c).
 *
 * Does NOT touch the BFS distance tables (unlike line_profile.c) — this
 * checks H13-I alone (a statement about toric classes, no d(pi) involved),
 * so it is not covered by the "no line-model BFS at n >= 10" guidance in
 * docs/notes/repo_state.md (that guidance is about d_walk/d_line via graph
 * search, not about I(pi)).
 *
 * Usage: h13i_exhaustive n [report_every_seconds]
 * Version h13i_exhaustive-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXN 13

static int I_of_pi(const int *pi, int n, int *rotbuf, int *posbuf) {
    /* rotbuf holds pi rotated by a: rotbuf[j] = pi[(a+j) % n] */
    memcpy(rotbuf, pi, n * sizeof(int));
    int inv00 = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (rotbuf[i] > rotbuf[j]) inv00++;
    int best = inv00;
    int inv_a = inv00;
    for (int a = 0; a < n; a++) {
        if (a > 0) {
            int front = rotbuf[0];
            memmove(rotbuf, rotbuf + 1, (n - 1) * sizeof(int));
            rotbuf[n - 1] = front;
            inv_a += (n - 1 - 2 * front);
        }
        if (inv_a < best) best = inv_a;
        for (int j = 0; j < n; j++) posbuf[rotbuf[j]] = j;
        int inv_b = inv_a;
        for (int b = 0; b < n - 1; b++) {
            int jstar = posbuf[b];
            inv_b += (n - 1 - 2 * jstar);
            if (inv_b < best) best = inv_b;
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
    if (argc < 2) { fprintf(stderr, "usage: %s n [report_every_seconds]\n", argv[0]); return 2; }
    int n = atoi(argv[1]);
    if (n < 2 || n > MAXN) { fprintf(stderr, "n out of range\n"); return 2; }
    int report_every = argc >= 3 ? atoi(argv[2]) : 30;

    int target = ((n - 1) * (n - 1)) / 4;
    int pi[MAXN], rotbuf[MAXN], posbuf[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;

    long total = 1;
    for (int k = 2; k <= n; k++) total *= k;

    long count = 0;
    int global_max = -1;
    long cnt_at_max = 0;
    long first_violation = -1;
    int violation_I = -1;
    int violation_pi[MAXN];

    time_t t0 = time(NULL), tlast = t0;
    do {
        int I = I_of_pi(pi, n, rotbuf, posbuf);
        if (I > global_max) { global_max = I; cnt_at_max = 1; }
        else if (I == global_max) cnt_at_max++;
        if (I > target && first_violation < 0) {
            first_violation = count;
            violation_I = I;
            memcpy(violation_pi, pi, n * sizeof(int));
        }
        count++;
        if (report_every > 0 && count % 2000000 == 0) {
            time_t tn = time(NULL);
            if (tn - tlast >= report_every) {
                fprintf(stderr, "n=%d progress %ld/%ld (%.1f%%) elapsed=%lds global_max=%d target=%d\n",
                        n, count, total, 100.0 * count / total, (long)(tn - t0), global_max, target);
                tlast = tn;
            }
        }
    } while (next_perm(pi, n));

    time_t t1 = time(NULL);
    printf("{\n");
    printf("  \"n\": %d,\n", n);
    printf("  \"total_permutations\": %ld,\n", total);
    printf("  \"checked\": %ld,\n", count);
    printf("  \"target_floor_(n-1)^2_over_4\": %d,\n", target);
    printf("  \"global_max_I\": %d,\n", global_max);
    printf("  \"count_at_max\": %ld,\n", cnt_at_max);
    printf("  \"H13_I_holds\": %s,\n", global_max <= target ? "true" : "false");
    if (first_violation >= 0) {
        printf("  \"first_violation_lexrank\": %ld,\n", first_violation);
        printf("  \"first_violation_I\": %d,\n", violation_I);
        printf("  \"first_violation_pi\": [");
        for (int i = 0; i < n; i++) printf("%d%s", violation_pi[i], i + 1 < n ? "," : "");
        printf("],\n");
    }
    printf("  \"elapsed_seconds\": %ld\n", (long)(t1 - t0));
    printf("}\n");
    return 0;
}
