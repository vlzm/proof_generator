/* toric_cut_scan.c -- exhaustive scan of the double-cut grid of every toric
 * class of permutations of Z_n.
 *
 * For pi in S_n and (a, b) in Z_n^2 let
 *      w^{a,b}_j = (pi(a + j) - b) mod n,   j = 0..n-1,
 *      F(a, b)   = inv(w^{a,b}).
 * The substitution (a, b) = (q + 1, q + 1 - c) is a bijection of Z_n^2, so this
 * is exactly the double-cut grid of docs/notes/h13_line_model.md and
 *      I(pi) = min_{a,b} F(a, b).
 * Every toric class contains a representative with pi(0) = 0, and I, mu and the
 * whole grid up to translation are class invariants, so the scan runs over the
 * (n-1)! representatives with pi(0) = 0.
 *
 * Reported per n:
 *   maxI      = max_pi I(pi)                      (H13-I: should be M = floor((n-1)^2/4))
 *   nI        = number of representatives with I(pi) = maxI
 *   argI      = lexicographically first such representative
 *   maxS      = max_pi sum_{b in Z_n} min_a F(a, b)   (C37: should be n*M)
 *   nS, argS  = count and lexicographically first argmax
 *   failI     = number of representatives with I(pi) > M
 *   failS     = number of representatives with sum_b min_a F(a,b) > n*M
 *
 * Usage: toric_cut_scan n [n2]      (scans n..n2, JSON lines on stdout)
 * Version toric_cut_scan-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXN 13

static int n;
static int perm[MAXN];
static int F[MAXN][MAXN];

/* F(0, b) by definition, then the rest by the increment recurrence
 *      F(a+1, b) = F(a, b) + n - 1 - 2 ((pi(a) - b) mod n),
 * which is Lemma 1 of docs/notes/h13i_verdict.md.  checks/check_C37.py
 * re-derives the same grid from the definition only. */
static void build_grid(void)
{
    for (int b = 0; b < n; b++) {
        int c = 0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++) {
                int vi = (perm[i] - b + n) % n, vj = (perm[j] - b + n) % n;
                if (vi > vj) c++;
            }
        F[0][b] = c;
    }
    for (int a = 0; a + 1 < n; a++)
        for (int b = 0; b < n; b++)
            F[a + 1][b] = F[a][b] + n - 1 - 2 * ((perm[a] - b + n) % n);
}

static void print_perm(const int *p)
{
    printf("[");
    for (int i = 0; i < n; i++) printf("%d%s", p[i], i + 1 < n ? "," : "");
    printf("]");
}

static void scan(void)
{
    int M = ((n - 1) * (n - 1)) / 4;
    long reps = 0, failI = 0, failS = 0;
    int maxI = -1; long nI = 0; int argI[MAXN];
    long maxS = -1, nS = 0; int argS[MAXN];
    clock_t t0 = clock();

    int a[MAXN], m = n - 1;
    for (int i = 0; i < m; i++) a[i] = i + 1;
    for (;;) {
        perm[0] = 0;
        for (int i = 0; i < m; i++) perm[i + 1] = a[i];
        build_grid();
        long s = 0; int I = 1 << 30;
        for (int b = 0; b < n; b++) {
            int mn = 1 << 30;
            for (int aa = 0; aa < n; aa++) if (F[aa][b] < mn) mn = F[aa][b];
            s += mn;
            if (mn < I) I = mn;
        }
        reps++;
        if (I > M) failI++;
        if (s > (long)n * M) failS++;
        if (I > maxI) { maxI = I; nI = 1; memcpy(argI, perm, sizeof perm); }
        else if (I == maxI) nI++;
        if (s > maxS) { maxS = s; nS = 1; memcpy(argS, perm, sizeof perm); }
        else if (s == maxS) nS++;

        int i = m - 2;
        while (i >= 0 && a[i] >= a[i + 1]) i--;
        if (i < 0) break;
        int j = m - 1;
        while (a[j] <= a[i]) j--;
        int t = a[i]; a[i] = a[j]; a[j] = t;
        for (int l = i + 1, r = m - 1; l < r; l++, r--) { t = a[l]; a[l] = a[r]; a[r] = t; }
    }
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("{\"n\": %d, \"M\": %d, \"nM\": %d, \"reps\": %ld, ", n, M, n * M, reps);
    printf("\"maxI\": %d, \"count_maxI\": %ld, \"arg_maxI\": ", maxI, nI);
    print_perm(argI);
    printf(", \"maxS\": %ld, \"count_maxS\": %ld, \"arg_maxS\": ", maxS, nS);
    print_perm(argS);
    printf(", \"fail_H13I\": %ld, \"fail_C37\": %ld, \"seconds\": %.2f}\n",
           failI, failS, secs);
    fflush(stdout);
}

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: %s n [n2]\n", argv[0]); return 2; }
    int n1 = atoi(argv[1]);
    int n2 = argc >= 3 ? atoi(argv[2]) : n1;
    if (n1 < 3 || n2 > MAXN - 1 || n2 < n1) return 2;
    for (n = n1; n <= n2; n++) scan();
    return 0;
}
