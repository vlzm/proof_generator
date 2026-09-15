/* toric_I_fast.c — exhaustive max_pi I(pi) over S_n, session 9 (H13-I).
 *
 * I(pi) = min over the n^2 double cuts (q, c) of inv(w), same operational
 * definition as experiments/line_profile.c: shift s = q+1, t = q+1-c (mod n),
 * w_j = (pi[(s+j) mod n] - t) mod n for j = 0..n-1, inv(w) = number of pairs
 * i<j with w_i>w_j.  Equivalently I(pi) = min over rotations of positions and
 * of values (the toric class of pi) of the linear inversion count.
 *
 * This tool replaces the O(n^2) per (q,c) x n^2 cuts = O(n^4)-per-permutation
 * scan of line_profile.c with an O(n^2)-per-permutation scan, using two exact
 * recurrences (verified against brute force for all pi, 2 <= n <= 7, and
 * against line_profile.c's max_pi I for 4 <= n <= 10; see
 * docs/notes/h13_line_model.md section 7 for the derivation and the proof
 * that these recurrences are exact identities, not heuristics):
 *
 *   inv(s+1, t) = inv(s, t) + (n-1) - 2*((pi[s]   - t) mod n)      (*)
 *   inv(s, t+1) = inv(s, t) + (n-1) - 2*((piinv[t] - s) mod n)     (**)
 *
 * (*) is the standard "move the front element to the back" bubble-sort
 * identity applied to the line at cut s: the element leaving the front has
 * value v = w_0(s,t) = (pi[s]-t) mod n; it was inverted with the v elements
 * smaller than it (all elements before it: none) and forms new inversions
 * with none of the larger ones ahead — concretely, moving it from front to
 * back changes inv by (n-1-v) - v = n-1-2v (loses v inversions with smaller
 * elements now behind it, gains n-1-v with the rest now ahead of it, net
 * counted directly on the two arrangements). (**) is the same identity
 * applied to the inverse line (rotating the value origin is, by inv(w) =
 * inv(w^{-1}), the same move on w^{-1}, whose "front element leaving" is the
 * position of the current minimum value, i.e. piinv[t] - s mod n).
 *
 * Algorithm per permutation: seed inv(0,0) directly (O(n^2), dominates the
 * O(n) per-row/per-column recurrence sweeps below only for the smallest n);
 * fill row t=0 for all s via (*); then for each s fill all t via (**).
 * Total O(n^2) per permutation, O(n^2 . n!) overall.
 *
 * Usage: toric_I_fast n
 * Prints: max_I, floor((n-1)^2/4), how many pi attain the max (and the rank
 * of the first one, Lehmer/lexicographic order), plus the full histogram of
 * I(pi) (sanity cross-check against experiments/line_profile.c at n <= 10).
 * Version toric_I_fast-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

static int inv0(const int *pi, int n) {
    int inv = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (pi[i] > pi[j]) inv++;
    return inv;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s n\n", argv[0]); return 2; }
    int n = atoi(argv[1]);
    if (n < 2 || n > MAXN) { fprintf(stderr, "n out of range [2, %d]\n", MAXN); return 2; }
    int pi[MAXN], piinv[MAXN], row0[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;
    long long total = 1;
    for (int k = 2; k <= n; k++) total *= k;
    int global_max = -1;
    long long argmax_rank = -1, cnt_at_max = 0;
    static long long hist[300];
    memset(hist, 0, sizeof hist);
    long long rank = 0;
    clock_t t0 = clock();
    do {
        for (int i = 0; i < n; i++) piinv[pi[i]] = i;
        int base = inv0(pi, n);
        row0[0] = base;
        for (int s = 0; s < n - 1; s++) {
            int v = pi[s]; /* (pi[s] - 0) mod n, t = 0 */
            row0[s + 1] = row0[s] + (n - 1) - 2 * v;
        }
        int best = base;
        for (int i = 1; i < n; i++) if (row0[i] < best) best = row0[i];
        for (int s = 0; s < n; s++) {
            int cur = row0[s];
            for (int t = 0; t < n - 1; t++) {
                int k = piinv[t] - s;
                if (k < 0) k += n;
                cur = cur + (n - 1) - 2 * k;
                if (cur < best) best = cur;
            }
        }
        hist[best]++;
        if (best > global_max) { global_max = best; argmax_rank = rank; cnt_at_max = 1; }
        else if (best == global_max) cnt_at_max++;
        rank++;
    } while (next_perm(pi, n));
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("{\"n\": %d, \"total\": %lld, \"max_I\": %d, \"floor_(n-1)^2/4\": %d, "
           "\"count_at_max\": %lld, \"argmax_rank\": %lld, \"time_s\": %.2f, \"histogram\": {",
           n, total, global_max, ((n - 1) * (n - 1)) / 4, cnt_at_max, argmax_rank, secs);
    int first = 1;
    for (int i = 0; i < 300; i++) if (hist[i]) {
        printf("%s\"%d\": %lld", first ? "" : ", ", i, hist[i]);
        first = 0;
    }
    printf("}}\n");
    return 0;
}
