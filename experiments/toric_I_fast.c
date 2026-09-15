/* toric_I_fast.c — exhaustive scan of I(pi) = min_{q,c} inv(w) (H13-I, see
 * docs/notes/h13_line_model.md) using an O(n^2)-per-permutation formula
 * instead of the O(n^4) brute force in experiments/line_profile.c. This is
 * what let session 9 push exhaustive verification of H13-I from n <= 10
 * (line_profile, session 8) to n <= 13.
 *
 * Double-cut convention (unchanged): w_j = pi(q+1+j) - (q+1-c) mod n,
 * j = 0..n-1. Writing s = q+1 mod n and t = (q+1-c) mod n (t ranges
 * independently over 0..n-1 as c does, for fixed s):
 *   F(s,t) := inv(w),  w_j = (pi((s+j) mod n) - t) mod n.
 *
 * Derivation (session 9, verified against brute force for all pi at
 * 2 <= n <= 7 and 2000 random pi at n = 8,9,10 in
 * experiments/verify_toric_I_fast.py / checks/check_C37.py):
 *
 *   F(0,0) = inv(pi)                                    (plain array inv).
 *   P(0) = F(0,0); P(s+1) = P(s) + (n-1-2*pi[s])         (rotate positions,
 *       drop front element pi[s] to the back: standard "rotate array by one"
 *       inversion-count delta, since values are a permutation of 0..n-1).
 *   For fixed s: F(s,0) = P(s);
 *     F(s,t+1) = F(s,t) + (n-1-2*sigma),  sigma = (rho(t)-s) mod n,
 *       rho = pi^{-1} (same rotate-array delta, applied to the VALUE
 *       sequence w(s,t), whose front element equals rho(t)-s taken mod n
 *       when t is the value about to become 0... concretely: incrementing
 *       t by 1 cyclically decrements every entry of w by 1 mod n, which
 *       moves the entry that was 0 (position rho(t)-s mod n in the row) to
 *       n-1; delta = n-1-2*(that position)).
 *   I(pi) = min_{s,t} F(s,t), computed by walking t = 0..n-1 for every s.
 *
 * This costs O(n) to build rho, O(n^2) (or O(n log n)) for inv(pi), O(n) for
 * all P(s), and O(n) per row (n rows) = O(n^2) total per permutation, vs.
 * O(n^4) for the brute-force n^2-cuts x O(n^2)-inversions approach.
 *
 * Usage: toric_I_fast n [progress]
 * Prints one JSON line: n, permutation count, floor((n-1)^2/4), max I found,
 * number of permutations exceeding the bound (should be 0 if H13-I holds),
 * the lexicographic-rank argmax and its permutation, and elapsed seconds.
 * With a second argument, progress lines go to stderr every 5,000,000 ranks.
 * Version toric_I_fast-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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
    if (argc < 2) { fprintf(stderr, "usage: %s n [progress]\n", argv[0]); return 2; }
    int n = atoi(argv[1]);
    if (n < 2 || n > MAXN) { fprintf(stderr, "n out of range [2, %d]\n", MAXN); return 2; }
    int bound = ((n - 1) * (n - 1)) / 4;

    long total = 1;
    for (int k = 2; k <= n; k++) total *= k;

    int pi[MAXN], rho[MAXN], P[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;

    int globalMax = -1;
    long argmaxRank = -1;
    int bestPi[MAXN];
    long rank = 0;
    long overCount = 0;
    time_t t0 = time(NULL);

    do {
        for (int i = 0; i < n; i++) rho[pi[i]] = i;
        int inv0 = 0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                if (pi[i] > pi[j]) inv0++;
        P[0] = inv0;
        for (int s = 0; s < n - 1; s++)
            P[s + 1] = P[s] + (n - 1 - 2 * pi[s]);

        int best = P[0];
        for (int s = 0; s < n; s++) {
            int F = P[s];
            if (F < best) best = F;
            for (int t = 0; t < n - 1; t++) {
                int sigma = rho[t] - s;
                if (sigma < 0) sigma += n;
                F += (n - 1 - 2 * sigma);
                if (F < best) best = F;
            }
        }

        if (best > globalMax) {
            globalMax = best;
            argmaxRank = rank;
            for (int i = 0; i < n; i++) bestPi[i] = pi[i];
        }
        if (best > bound) overCount++;

        rank++;
        if (argc >= 3 && rank % 5000000 == 0) {
            fprintf(stderr, "n=%d rank=%ld/%ld (%.1f%%) elapsed=%lds maxI=%d\n",
                    n, rank, total, 100.0 * rank / total, (long)(time(NULL) - t0), globalMax);
        }
    } while (next_perm(pi, n));

    printf("{\"n\": %d, \"count\": %ld, \"bound_floor_(n-1)^2/4\": %d, \"max_I\": %d, "
           "\"exceeds_count\": %ld, \"argmax_rank\": %ld, \"argmax_pi\": [",
           n, total, bound, globalMax, overCount, argmaxRank);
    for (int i = 0; i < n; i++) printf("%d%s", bestPi[i], i + 1 < n ? "," : "");
    printf("], \"elapsed_s\": %ld}\n", (long)(time(NULL) - t0));
    return 0;
}
