/* h13i_avgq.c — exhaustive check of a candidate route to H13-I
 * (`I(pi) = min_{q,c} inv(w) <= floor((n-1)^2/4)`, see docs/notes/h13_line_model.md §6).
 *
 * For a permutation pi of {0,...,n-1} and a position-cut q, let v be pi
 * rotated to start right after q (v_j = pi[(q+1+j) mod n]) and let
 * M(q) = min_c inv(shift_c(v)) (shift_c(v)_j = (v_j - c) mod n). M(q) is
 * computed in O(n^2) via the "cycle walk": inv(v) at c=0, then a single
 * pass c = 0..n-1 where each step changes inv by (n-1-2*p_c), p_c the
 * position of value c in v (standard fact for cyclic value-relabelling).
 *
 * Two modes:
 *   avgq   (default) — for every pi, compute avg_q M(q) and compare to
 *          floor((n-1)^2/4). Tests the conjecture "H13-I-avg":
 *          avg_q M(q) <= floor((n-1)^2/4) for all pi, which would prove
 *          H13-I by averaging (min_q M(q) <= avg_q M(q)).
 *   singleq — for every pi, compute M(0) alone (q fixed at 0, i.e. only the
 *          value-shift c is optimized) and compare to floor((n-1)^2/4).
 *          Rules out "value-shift alone suffices" as a proof route.
 *
 * Usage: h13i_avgq n [avgq|singleq]
 * Version h13i_avgq-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 12

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

/* M(v) = min_c inv(shift_c(v)) for a permutation array v of size n. */
static int min_c_inv(const int *v, int n) {
    static int p[MAXN];
    for (int j = 0; j < n; j++) p[v[j]] = j;
    int inv0 = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (v[i] > v[j]) inv0++;
    int S = 0, minS = 0;
    for (int c = 0; c < n; c++) {
        S += (n - 1) - 2 * p[c];
        if (S < minS) minS = S;
    }
    return inv0 + minS;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s n [avgq|singleq]\n", argv[0]); return 2; }
    int n = atoi(argv[1]);
    if (n < 2 || n > MAXN) return 2;
    const char *mode = argc >= 3 ? argv[2] : "avgq";
    int bound = (n - 1) * (n - 1) / 4;

    int pi[MAXN], v[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;

    long best_num = -1, argmax_rank = -1, rank = 0, violations = 0, total_num = 0;

    do {
        long value; /* n * M for avgq (so we can report exact rationals), or M for singleq */
        if (strcmp(mode, "singleq") == 0) {
            value = min_c_inv(pi, n);
            total_num += value;
            if (value > best_num) { best_num = value; argmax_rank = rank; }
            if (value > bound) violations++;
        } else {
            long sum = 0;
            for (int q = 0; q < n; q++) {
                for (int j = 0; j < n; j++) v[j] = pi[(q + 1 + j) % n];
                sum += min_c_inv(v, n);
            }
            total_num += sum;
            if (sum > best_num) { best_num = sum; argmax_rank = rank; }
            if (sum > (long)bound * n) violations++;
        }
        rank++;
    } while (next_perm(pi, n));

    if (strcmp(mode, "singleq") == 0) {
        printf("n=%d mode=singleq bound=%d max_M=%ld argmax_rank=%ld violations=%ld total_perms=%ld avg_M=%.6f\n",
               n, bound, best_num, argmax_rank, violations, rank, (double)total_num / rank);
    } else {
        printf("n=%d mode=avgq bound=%d max_avg_q_M=%.6f (=%ld/%d) argmax_rank=%ld violations=%ld total_perms=%ld\n",
               n, bound, (double)best_num / n, best_num, n, argmax_rank, violations, rank);
    }
    return 0;
}
