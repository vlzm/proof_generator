/* h13i_avgq-1.0
 *
 * H13-I candidate reformulation. For pi in S_n and a position-rotation q,
 * define the q-rotated array b^(q)_j = pi((q+j) mod n), j = 0..n-1, and
 *   m(q) = min_{c in Z_n} inv( (b^(q)_j - c) mod n : j = 0..n-1 ).
 * I(pi) = min_q m(q) is the minimum inversions over all n^2 double cuts
 * (position rotation q, value rotation c); H13-I claims I(pi) <=
 * floor((n-1)^2/4) for every pi.
 *
 * This tool tests the stronger AVERAGING claim
 *   sum_{q=0}^{n-1} m(q) <= n * floor((n-1)^2/4)               (*)
 * which implies H13-I by pigeonhole (some q has m(q) <= the average).
 * Unlike averaging inv(q,c) directly over all n^2 cuts (refuted: the
 * average exceeds the bound already on sigma_n, id), here c is minimized
 * out first for each q, and only then q is averaged.
 *
 * inv_q(0) is computed in O(n^2); the rest of the row inv_q(1..n-1) is
 * obtained in O(n) via the exact recursion
 *   inv_q(c+1) - inv_q(c) = n - 1 - 2*pos_q(c),  pos_q = inverse of b^(q)
 * (moving the element valued c, at rank 0, to rank n-1 changes its own
 * contribution from pos_q(c) inversions to n-1-pos_q(c); every other
 * element's relative order is unchanged). This gives O(n) per q after an
 * O(n^2) start, O(n^3) per permutation, O(n! n^3) for exhaustive n.
 *
 * Modes:
 *   exhaustive N   : all N! permutations, N <= ~10
 *   sample N T     : identity, reversal + rotations, affine (a,b), T
 *                    random permutations, report worst sum_q m(q)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int n;

static long inv_count(const int *seq, int len) {
    long c = 0;
    for (int i = 0; i < len; i++)
        for (int j = i + 1; j < len; j++)
            if (seq[i] > seq[j]) c++;
    return c;
}

static long m_of_q(const int *pi_, int q) {
    int rotated[4096], posq[4096];
    for (int j = 0; j < n; j++) rotated[j] = pi_[(q + j) % n];
    for (int j = 0; j < n; j++) posq[rotated[j]] = j;
    long cur = inv_count(rotated, n);
    long best = cur;
    for (int c = 0; c < n - 1; c++) {
        cur += (n - 1 - 2 * posq[c]);
        if (cur < best) best = cur;
    }
    return best;
}

static long sum_over_q(const int *pi_) {
    long s = 0;
    for (int q = 0; q < n; q++) s += m_of_q(pi_, q);
    return s;
}

static long bound_n(void) { return (long)n * ((long)(n - 1) * (n - 1) / 4); }

static long worst_sum = -1;
static int worst_perm[4096];
static long checked = 0;

static void consider(const int *pi_) {
    long s = sum_over_q(pi_);
    if (s > worst_sum) { worst_sum = s; memcpy(worst_perm, pi_, sizeof(int) * n); }
    checked++;
}

static int cur[4096];
static void permute(int k) {
    if (k == n) { consider(cur); return; }
    for (int i = k; i < n; i++) {
        int t = cur[k]; cur[k] = cur[i]; cur[i] = t;
        permute(k + 1);
        t = cur[k]; cur[k] = cur[i]; cur[i] = t;
    }
}

static void run_exhaustive(void) {
    for (int i = 0; i < n; i++) cur[i] = i;
    permute(0);
    long bn = bound_n();
    printf("exhaustive n=%d checked=%ld bound_n=%ld worst_sum=%ld exceeds=%s avg=%.4f target_per_q=%ld worst_perm=",
           n, checked, bn, worst_sum, worst_sum > bn ? "YES" : "no",
           (double)worst_sum / n, (long)(n - 1) * (n - 1) / 4);
    for (int i = 0; i < n; i++) printf("%d%s", worst_perm[i], i + 1 < n ? "," : "\n");
}

static void run_sample(int trials) {
    int *pi_ = malloc(sizeof(int) * n);
    long bn = bound_n();
    long global_worst = -1;
    int any_exceed = 0;

    for (int i = 0; i < n; i++) pi_[i] = i;
    { long s = sum_over_q(pi_); if (s > global_worst) global_worst = s; }

    for (int i = 0; i < n; i++) pi_[i] = n - 1 - i;
    { long s = sum_over_q(pi_); if (s > global_worst) global_worst = s; if (s > bn) any_exceed = 1; }

    for (int r = 0; r < n; r += (n / 11 + 1)) {
        for (int i = 0; i < n; i++) pi_[i] = ((n - 1 - i) + r) % n;
        long s = sum_over_q(pi_);
        if (s > global_worst) global_worst = s;
        if (s > bn) any_exceed = 1;
    }

    for (int a = 2; a < n && a < 14; a++) {
        int g0 = a, g1 = n, t;
        while (g1) { t = g0 % g1; g0 = g1; g1 = t; }
        if (g0 != 1) continue;
        for (int b = 0; b < n; b += (n / 5 + 1)) {
            for (int i = 0; i < n; i++) pi_[i] = ((long)a * i + b) % n;
            long s = sum_over_q(pi_);
            if (s > global_worst) global_worst = s;
            if (s > bn) any_exceed = 1;
        }
    }

    srand(20260913);
    for (int trial = 0; trial < trials; trial++) {
        for (int i = 0; i < n; i++) pi_[i] = i;
        for (int i = n - 1; i > 0; i--) { int j = rand() % (i + 1); int t = pi_[i]; pi_[i] = pi_[j]; pi_[j] = t; }
        long s = sum_over_q(pi_);
        if (s > global_worst) global_worst = s;
        if (s > bn) any_exceed = 1;
    }

    printf("sample n=%d trials=%d bound_n=%ld worst_sum=%ld exceeds=%s avg=%.4f target_per_q=%ld\n",
           n, trials, bn, global_worst, any_exceed ? "YES" : "no",
           (double)global_worst / n, (long)(n - 1) * (n - 1) / 4);
    free(pi_);
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s exhaustive N | %s sample N trials\n", argv[0], argv[0]);
        return 1;
    }
    n = atoi(argv[2]);
    if (strcmp(argv[1], "exhaustive") == 0) {
        run_exhaustive();
    } else if (strcmp(argv[1], "sample") == 0) {
        int trials = argc > 3 ? atoi(argv[3]) : 200;
        run_sample(trials);
    } else {
        fprintf(stderr, "unknown mode %s\n", argv[1]);
        return 1;
    }
    return 0;
}
