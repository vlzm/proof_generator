/*
 * Exhaustive scan of I(pi) = min over double cuts (q, c) of inv(w) for all
 * permutations pi of size n (H13-I / C33, docs/notes/h13_line_model.md).
 *
 * A double cut is a choice of position-rotation a and value-rotation b
 * (n^2 choices); w_j = (pi((a+j) mod n) - b) mod n, j = 0..n-1, and
 * inv(w) counts pairs j < k with w_j > w_k. I(pi) = min_{a,b} inv(w).
 *
 * Naive evaluation costs O(n^2) cuts * O(n^2) inversions = O(n^4) per pi,
 * making exhaustive search infeasible beyond n ~ 9-10. This scanner uses
 * the closed-form gradient identity (C37, checks/check_C37.py):
 *
 *   g(a+1, b) - g(a, b) = n - 1 - 2 * ((pi(a) - b) mod n)
 *   g(a, b+1) - g(a, b) = n - 1 - 2 * ((pi^{-1}(b) - a) mod n)
 *
 * where g(a, b) = inv(w) for cut (a, b). This lets the full n x n grid of
 * g values be filled in O(n^2) per permutation (one base inv(pi) computed
 * directly, then O(1) per grid cell), giving O(n! * n^2) total instead of
 * O(n! * n^4). Verified: worst-case (max over pi of I(pi)) equals the
 * conjectured bound floor((n-1)^2/4) for 4 <= n <= 12 (this run extends
 * the prior exhaustive range 4 <= n <= 10 to 4 <= n <= 12).
 *
 * Usage: ./toric_inversions_scan N
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int n;
static int pi_[16], invpi[16];
static int g[16][16];
static long long worst;
static int worst_perm[16];

static int inv_count(const int *a, int n) {
    int c = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (a[i] > a[j]) c++;
    return c;
}

static void process(void) {
    for (int i = 0; i < n; i++) invpi[pi_[i]] = i;

    g[0][0] = inv_count(pi_, n);
    for (int a = 0; a < n - 1; a++) {
        int step = n - 1 - 2 * pi_[a];
        g[a + 1][0] = g[a][0] + step;
    }
    for (int a = 0; a < n; a++) {
        for (int b = 0; b < n - 1; b++) {
            int jstar = invpi[b] - a;
            jstar %= n;
            if (jstar < 0) jstar += n;
            int step = n - 1 - 2 * jstar;
            g[a][b + 1] = g[a][b] + step;
        }
    }

    int mn = g[0][0];
    for (int a = 0; a < n; a++)
        for (int b = 0; b < n; b++)
            if (g[a][b] < mn) mn = g[a][b];

    if (mn > worst) {
        worst = mn;
        memcpy(worst_perm, pi_, sizeof(int) * n);
    }
}

static void permute(int k) {
    if (k == n) { process(); return; }
    for (int i = k; i < n; i++) {
        int t = pi_[k]; pi_[k] = pi_[i]; pi_[i] = t;
        permute(k + 1);
        t = pi_[k]; pi_[k] = pi_[i]; pi_[i] = t;
    }
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s N\n", argv[0]); return 1; }
    n = atoi(argv[1]);
    if (n < 1 || n > 13) { fprintf(stderr, "N out of supported range (1..13)\n"); return 1; }
    for (int i = 0; i < n; i++) pi_[i] = i;
    worst = -1;
    permute(0);
    long long bound = (long long)(n - 1) * (n - 1) / 4;
    printf("n=%d bound=%lld worst=%lld %s\n", n, bound, worst, worst <= bound ? "OK" : "FAIL");
    printf("worst_perm:");
    for (int i = 0; i < n; i++) printf(" %d", worst_perm[i]);
    printf("\n");
    return 0;
}
