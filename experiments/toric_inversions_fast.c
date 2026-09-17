/* toric_inversions_fast.c — I(pi) = min over double cuts (q, c) of the number
 * of inversions of the relabelled line (same definition as
 * experiments/line_profile.c and docs/notes/h13_line_model.md), for every
 * permutation of n (Lehmer / lexicographic order).
 *
 * Speedup over line_profile.c: instead of recomputing inversions from
 * scratch for each of the n^2 cuts (O(n^2) per cut, O(n^4) per permutation),
 * this uses the exact recursions
 *   inv(a+1, b) = inv(a, b) + (n-1) - 2*((pi[a] - b) mod n)
 *   inv(0, b+1) = inv(0, b) + (n-1) - 2*((pinv[b] - 0) mod n)
 * (moving the head element of the a-rotation, resp. the value-b element, to
 * the back of its cycle) to fill the whole n x n table of inv(a, b) in O(n^2)
 * time after one O(n^2) base computation of inv(0, 0). Total O(n! * n^2).
 *
 * Purpose: extend the exhaustive check of H13-I (max_pi I(pi) =
 * floor((n-1)^2/4), attained only on the n reflections) from 4 <= n <= 10
 * (line_profile.c, session 8) to 4 <= n <= 12 within the session-9 budget.
 * Cross-checked against line_profile.c's reported max_I and reflection count
 * at 4 <= n <= 10 (identical), and against a brute-force O(n^4) Python
 * re-implementation at 4 <= n <= 8 (identical, including the exact set of
 * extremal permutations = the n reflections).
 *
 * Usage: toric_inversions_fast n
 * Output: one line "n=.. total=.. maxI=.. bound=.. countMax=.."
 * Version toric_inversions_fast-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 14

static int n;
static int pi_[MAXN], pinv[MAXN];

static int inv0(void) {
    int c = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (pi_[i] > pi_[j]) c++;
    return c;
}

static int next_perm(int *a, int len) {
    int i = len - 2;
    while (i >= 0 && a[i] >= a[i + 1]) i--;
    if (i < 0) return 0;
    int j = len - 1;
    while (a[j] <= a[i]) j--;
    int t = a[i]; a[i] = a[j]; a[j] = t;
    for (int l = i + 1, r = len - 1; l < r; l++, r--) { t = a[l]; a[l] = a[r]; a[r] = t; }
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s n\n", argv[0]); return 2; }
    n = atoi(argv[1]);
    if (n < 2 || n > MAXN) { fprintf(stderr, "n out of range\n"); return 2; }
    for (int i = 0; i < n; i++) pi_[i] = i;
    long total = 1;
    for (int k = 2; k <= n; k++) total *= k;

    int maxI = -1;
    long countMax = 0;
    int row0[MAXN];
    do {
        for (int i = 0; i < n; i++) pinv[pi_[i]] = i;
        int base = inv0();
        row0[0] = base;
        for (int b = 0; b < n - 1; b++) {
            int r = pinv[b];
            row0[b + 1] = row0[b] + (n - 1) - 2 * r;
        }
        int best = 1000000;
        for (int b = 0; b < n; b++) {
            int cur = row0[b];
            if (cur < best) best = cur;
            for (int a = 0; a < n - 1; a++) {
                int r = pi_[a] - b;
                if (r < 0) r += n;
                cur = cur + (n - 1) - 2 * r;
                if (cur < best) best = cur;
            }
        }
        if (best > maxI) { maxI = best; countMax = 1; }
        else if (best == maxI) { countMax++; }
    } while (next_perm(pi_, n));

    printf("n=%d total=%ld maxI=%d bound=%d countMax=%ld\n",
           n, total, maxI, ((n - 1) * (n - 1)) / 4, countMax);
    return 0;
}
