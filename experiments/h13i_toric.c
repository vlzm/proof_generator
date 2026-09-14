/* h13i_toric.c -- exhaustive max_pi I(pi) for hypothesis H13-I (PLAN.md H13,
 * docs/notes/h13_line_model.md), extending the range beyond experiments/
 * line_profile.c (which needs the distance table and is capped by the table
 * size). H13-I only concerns I(pi), so no distance table is needed here.
 *
 * I(pi) = min over position-cut q in Z_n and value-shift c in Z_n of the
 * number of inversions of the relabelled line
 *   w_j = (pi[(q+1+j) mod n] - shift) mod n,  shift = q+1-c,  j = 0..n-1
 * (same definition as line_profile.c / line_model.py, verified to agree with
 * the brute-force O(n^2) x O(n) enumeration on random permutations up to
 * n = 8, see data/runs/h13i_search/).
 *
 * For fixed q, instead of looping over all n shifts and recomputing
 * inversions from scratch (O(n^3) per q), use: for a pair of line positions
 * j < k with values a_j, a_k (a = pi rotated to start at q+1), the set of
 * shifts s causing an inversion is a contiguous circular arc of length n-d
 * (d = (a_j - a_k) mod n) ending at a_k; accumulate all C(n,2) pairs with a
 * circular difference array in O(n) extra space, then prefix-sum once to
 * read off inv(s) for every s. Total per pi: O(n^3) (O(n^2) pairs per q,
 * n values of q); no BFS, no distance table.
 *
 * Usage: h13i_toric n [skip_mod skip_res]
 *   n            permutation size, 2 <= n <= MAXN
 *   skip_mod,res  optional: only process permutations with
 *                 (lexicographic_rank - 1) % skip_mod == skip_res, for
 *                 sharding n! across processes.
 * Prints one line: n, permutations checked, target = floor((n-1)^2/4),
 * observed max_pi I(pi), the lexicographic rank of one argmax, the argmax
 * itself, and ok/COUNTEREXAMPLE.
 * Version h13i_toric-1.0 (session 9, 14.09.2026).
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
    if (argc < 2) { fprintf(stderr, "usage: %s n [skip_mod skip_res]\n", argv[0]); return 2; }
    int n = atoi(argv[1]);
    if (n < 2 || n > MAXN) { fprintf(stderr, "n out of range 2..%d\n", MAXN); return 2; }
    long skip_mod = argc >= 4 ? atol(argv[2]) : 1;
    long skip_res = argc >= 4 ? atol(argv[3]) : 0;

    int pi[MAXN], a[MAXN];
    int diff[MAXN + 1];
    for (int i = 0; i < n; i++) pi[i] = i;

    int global_best = -1;
    long argbest_rank = -1;
    int argbest[MAXN];
    long total = 0, rank = 0;

    do {
        rank++;
        if (skip_mod != 1 && (rank - 1) % skip_mod != skip_res) continue;
        total++;
        int per_pi_best = n * n;
        for (int q = 0; q < n; q++) {
            for (int j = 0; j < n; j++) a[j] = pi[(q + 1 + j) % n];
            memset(diff, 0, sizeof(int) * (n + 1));
            for (int j = 0; j < n; j++) {
                int aj = a[j];
                for (int k = j + 1; k < n; k++) {
                    int ak = a[k];
                    int d = aj - ak; if (d < 0) d += n;
                    int length = n - d;
                    int start = ak - length + 1;
                    start %= n; if (start < 0) start += n;
                    int end = start + length;
                    if (end <= n) {
                        diff[start] += 1;
                        diff[end] -= 1;
                    } else {
                        diff[start] += 1;
                        diff[n] -= 1;
                        diff[0] += 1;
                        diff[end - n] -= 1;
                    }
                }
            }
            int cur = 0;
            for (int s = 0; s < n; s++) {
                cur += diff[s];
                if (cur < per_pi_best) per_pi_best = cur;
            }
        }
        if (per_pi_best > global_best) {
            global_best = per_pi_best;
            argbest_rank = rank;
            memcpy(argbest, pi, sizeof(int) * n);
        }
    } while (next_perm(pi, n));

    int target = (n - 1) * (n - 1) / 4;
    printf("n=%d checked=%ld target=%d max_I=%d rank_of_argmax=%ld argmax=(", n, total, target, global_best, argbest_rank);
    for (int i = 0; i < n; i++) printf("%d%s", argbest[i], i + 1 < n ? "," : "");
    printf(") %s\n", global_best > target ? "COUNTEREXAMPLE" : "ok");
    return 0;
}
