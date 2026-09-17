/* line_toric_max_i.c -- exhaustive max_pi I(pi) over S_n (H13-I, session 9).
 *
 * I(pi) = min_{q,c} inv(w) where w_j = pi[(q+1+j) mod n] - (q+1-c) mod n,
 * j = 0..n-1 (same double cut as experiments/line_model.py and
 * experiments/line_profile.c).  This file is an INDEPENDENT reimplementation
 * (fresh code, not sharing logic with line_profile.c) that computes only
 * max_pi I(pi) -- not the joint profile against d(pi) -- using an O(n)
 * incremental update over the value-shift c, instead of recomputing
 * inv(w) from scratch (O(n^2)) for every one of the n choices of c.
 *
 * Incremental step (fixed q, hence fixed rotated line v[j] = pi[(q+1+j) mod n]):
 * let y_j(shift) = (v[j] - shift) mod n.  Going shift -> shift+1, every y_j
 * decreases by 1 except the one position p with v[p] == shift, whose value
 * wraps from 0 to n-1.  A short case check shows
 *     inv(shift+1) - inv(shift) = n - 1 - 2*p
 * (the wrapped element loses its p inversions with the p elements now to its
 * left and gains n-1-p new inversions with the elements to its right).  This
 * turns the per-permutation cost from O(n) choices of q times O(n) choices of
 * c times O(n^2) inversions = O(n^4) into O(n) * (O(n^2) once at shift=0 +
 * O(n) for the remaining n-1 shifts) = O(n^3).
 *
 * Validated (checks/check_H13I.py) against a fresh, independent O(n^4) Python
 * brute force for 4 <= n <= 8, and against the existing exhaustive results in
 * data/runs/line_profile/profile_n{4..10}.json (max_I column) for 4 <= n <= 10.
 *
 * Usage: line_toric_max_i n [report_every]
 * Prints one JSON object: n, count (=n!), max_I, floor((n-1)^2/4), the rank
 * (lexicographic order of next_perm) and value of one permutation achieving
 * max_I, and the full histogram of I over all pi.
 * Version line_toric_max_i-1.0.
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
    if (argc < 2) { fprintf(stderr, "usage: %s n [report_every]\n", argv[0]); return 2; }
    int n = atoi(argv[1]);
    long report_every = (argc >= 3) ? atol(argv[2]) : 0L;
    if (n < 2 || n > MAXN) { fprintf(stderr, "n out of range (2..%d)\n", MAXN); return 2; }

    int pi[MAXN], v[MAXN], pos[MAXN], y[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;

    long total = 1;
    for (int k = 2; k <= n; k++) total *= k;

    int global_best = -1;
    long global_best_rank = -1;
    int global_best_perm[MAXN];
    long cnt_by_I[200];
    memset(cnt_by_I, 0, sizeof cnt_by_I);

    long rank = 0;
    do {
        int pi_best = 1000000;
        for (int q = 0; q < n; q++) {
            for (int j = 0; j < n; j++) v[j] = pi[(q + 1 + j) % n];
            for (int j = 0; j < n; j++) pos[v[j]] = j;
            int inv = 0;
            for (int j = 0; j < n; j++) y[j] = v[j];
            for (int i2 = 0; i2 < n; i2++)
                for (int j2 = i2 + 1; j2 < n; j2++)
                    if (y[i2] > y[j2]) inv++;
            if (inv < pi_best) pi_best = inv;
            for (int shift = 1; shift < n; shift++) {
                int p = pos[shift - 1];
                inv += (n - 1 - 2 * p);
                if (inv < pi_best) pi_best = inv;
            }
        }
        if (pi_best >= 0 && pi_best < 200) cnt_by_I[pi_best]++;
        if (pi_best > global_best) {
            global_best = pi_best;
            global_best_rank = rank;
            memcpy(global_best_perm, pi, sizeof(int) * n);
        }
        rank++;
        if (report_every > 0 && rank % report_every == 0)
            fprintf(stderr, "n=%d progress %ld/%ld best_so_far=%d\n", n, rank, total, global_best);
    } while (next_perm(pi, n));

    printf("{\"n\": %d, \"count\": %ld, \"max_I\": %d, \"floor_(n-1)2_4\": %d, "
           "\"argmax_rank\": %ld, \"argmax\": [",
           n, total, global_best, (n - 1) * (n - 1) / 4, global_best_rank);
    for (int i = 0; i < n; i++) printf("%s%d", i ? "," : "", global_best_perm[i]);
    printf("], \"cnt_by_I\": {");
    int first = 1;
    for (int i = 0; i < 200; i++) if (cnt_by_I[i]) {
        printf("%s\"%d\": %ld", first ? "" : ", ", i, cnt_by_I[i]);
        first = 0;
    }
    printf("}}\n");
    return 0;
}
