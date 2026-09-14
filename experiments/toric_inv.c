/* toric_inv.c — statistics of the toric class of a permutation of Z_n.
 *
 * For a permutation pi of Z_n and a double cut (u, v) the line is
 *     w_j = (pi(u + j) - v) mod n,   j = 0..n-1,
 * and inv(u, v) is its number of linear inversions.  Quantities:
 *     I(pi)     = min_{u,v} inv(u,v)                     (H13-I: I <= floor((n-1)^2/4))
 *     rowmin(u) = min_v inv(u,v)
 *     RA(pi)    = sum_u rowmin(u)                        (row-average form: RA <= n*floor((n-1)^2/4))
 *     Ibar(pi)  = (1/n^2) sum_{u,v} inv(u,v)             (reported as n^2 * Ibar, an integer)
 *
 * Modes:
 *   all n            exhaustive over one representative per toric class
 *                    (all pi with pi(0) = 0; every class is met exactly n times)
 *   affine nmin nmax all affine pi(i) = a*i + b with gcd(a,n) = 1, all a and all b
 *                    (b only shifts the value cut, so it changes nothing; checked anyway)
 *   rand n trials seed   random permutations (Fisher-Yates, xorshift64*)
 *   climb n trials seed  hill climbing on I(pi) by random transpositions
 *
 * Output: JSON lines on stdout.  Version toric_inv-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXN 200

static int n;
static int inv_row[MAXN * MAXN]; /* inv(u,v) as a flat n x n table */

static int n_good_rows; /* rows u with min_v inv(u,v) <= floor((n-1)^2/4) */

/* fill inv[] for pi; returns I = min, *psum_rowmin = sum_u min_v inv(u,v),
 * *pIbar_num = sum_{u,v} inv(u,v) */
static int toric(const int *pi, const int *sigma, long *psum_rowmin, long *pIbar_num,
                 int *pbest_u, int *pbest_v) {
    int inv0 = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (pi[i] > pi[j]) inv0++;
    /* row u = 0 across v: inv(0,v+1) = inv(0,v) + n-1-2*sigma(v) */
    inv_row[0] = inv0;
    for (int v = 0; v + 1 < n; v++)
        inv_row[v + 1] = inv_row[v] + n - 1 - 2 * sigma[v];
    /* down the columns: inv(u+1,v) = inv(u,v) + n-1-2*((pi(u)-v) mod n) */
    for (int u = 0; u + 1 < n; u++) {
        const int *src = inv_row + (long)u * n;
        int *dst = inv_row + (long)(u + 1) * n;
        int p = pi[u];
        for (int v = 0; v < n; v++) {
            int a = p - v; if (a < 0) a += n;
            dst[v] = src[v] + n - 1 - 2 * a;
        }
    }
    int I = 1 << 30; long sum_rowmin = 0, tot = 0;
    int bu = 0, bv = 0, good = 0, bnd = ((n - 1) * (n - 1)) / 4;
    for (int u = 0; u < n; u++) {
        const int *r = inv_row + (long)u * n;
        int rm = 1 << 30, rv = 0;
        for (int v = 0; v < n; v++) { tot += r[v]; if (r[v] < rm) { rm = r[v]; rv = v; } }
        sum_rowmin += rm;
        if (rm <= bnd) good++;
        if (rm < I) { I = rm; bu = u; bv = rv; }
    }
    n_good_rows = good;
    *psum_rowmin = sum_rowmin; *pIbar_num = tot; *pbest_u = bu; *pbest_v = bv;
    return I;
}

static void inverse(const int *pi, int *sigma) {
    for (int i = 0; i < n; i++) sigma[pi[i]] = i;
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

static void print_perm(const int *pi) {
    printf("[");
    for (int i = 0; i < n; i++) printf("%d%s", pi[i], i + 1 < n ? "," : "");
    printf("]");
}

static int bound(int nn) { return ((nn - 1) * (nn - 1)) / 4; }

static unsigned long long rng_state;
static unsigned long long rnd(void) {
    unsigned long long x = rng_state;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    rng_state = x;
    return x * 2685821657736338717ULL;
}

static void report(const char *mode, long long checked, int maxI, const int *argmaxI,
                   long maxRA, const int *argmaxRA, long nviol, const int *viol,
                   double secs) {
    printf("{\"mode\":\"%s\",\"n\":%d,\"bound\":%d,\"checked\":%lld,"
           "\"maxI\":%d,\"argmaxI\":", mode, n, bound(n), checked, maxI);
    print_perm(argmaxI);
    printf(",\"maxRA\":%ld,\"n_times_bound\":%d,\"argmaxRA\":", maxRA, n * bound(n));
    print_perm(argmaxRA);
    printf(",\"violations_I\":%ld", nviol);
    if (nviol) { printf(",\"first_violation\":"); print_perm(viol); }
    printf(",\"seconds\":%.1f}\n", secs);
    fflush(stdout);
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s all n | affine nmin nmax | rand n trials seed | climb n trials seed\n", argv[0]);
        return 2;
    }
    const char *mode = argv[1];
    int pi[MAXN], sigma[MAXN], argI[MAXN], argRA[MAXN], viol[MAXN];
    long sum_rowmin, tot; int bu, bv;
    clock_t t0 = clock();

    if (!strcmp(mode, "all")) {
        n = atoi(argv[2]);
        if (n < 2 || n > 13) return 2;   /* 13 takes ~3.5 min; 14 would take ~50 min */
        int maxI = -1; long maxRA = -1, nviol = 0, nviolRA = 0; long long checked = 0;
        long eqI = 0, eqRA = 0; int minGood = 1 << 30; int argGood[MAXN];
        long maxTot = -1; int argTot[MAXN]; memset(argTot, 0, sizeof argTot);
        memset(argI, 0, sizeof argI); memset(argRA, 0, sizeof argRA); memset(viol, 0, sizeof viol);
        memset(argGood, 0, sizeof argGood);
        int tail[MAXN];
        for (int i = 0; i + 1 < n; i++) tail[i] = i + 1;
        do {
            pi[0] = 0;
            for (int i = 0; i + 1 < n; i++) pi[i + 1] = tail[i];
            inverse(pi, sigma);
            int I = toric(pi, sigma, &sum_rowmin, &tot, &bu, &bv);
            checked++;
            if (I > maxI) { maxI = I; memcpy(argI, pi, sizeof pi); }
            if (sum_rowmin > maxRA) { maxRA = sum_rowmin; memcpy(argRA, pi, sizeof pi); }
            if (tot > maxTot) { maxTot = tot; memcpy(argTot, pi, sizeof pi); }
            if (I == bound(n)) eqI++;
            if (sum_rowmin == (long)n * bound(n)) eqRA++;
            if (n_good_rows < minGood) { minGood = n_good_rows; memcpy(argGood, pi, sizeof pi); }
            if (I > bound(n)) { if (!nviol) memcpy(viol, pi, sizeof pi); nviol++; }
            if (sum_rowmin > (long)n * bound(n)) nviolRA++;
        } while (next_perm(tail, n - 1));
        double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
        report("all", checked, maxI, argI, maxRA, argRA, nviol, viol, secs);
        printf("{\"mode\":\"all\",\"n\":%d,\"violations_RA\":%ld,\"equality_I\":%ld,"
               "\"equality_RA\":%ld,\"min_good_rows\":%d,\"argmin_good_rows\":", n, nviolRA, eqI, eqRA, minGood);
        print_perm(argGood);
        printf(",\"max_n2_Ibar\":%ld,\"reflection_n2_Ibar\":%ld,\"argmax_Ibar\":", maxTot,
               (long)n * n * (n - 1) * (2 * n - 1) / 6);
        print_perm(argTot);
        printf("}\n");
    } else if (!strcmp(mode, "affine")) {
        int nmin = atoi(argv[2]), nmax = atoi(argv[3]);
        for (n = nmin; n <= nmax; n++) {
            int maxI = -1; long maxRA = -1, nviol = 0, nviolRA = 0; long long checked = 0;
            memset(argI, 0, sizeof argI); memset(argRA, 0, sizeof argRA); memset(viol, 0, sizeof viol);
            for (int a = 1; a < n; a++) {
                int g = a, h = n; while (h) { int t = g % h; g = h; h = t; }
                if (g != 1) continue;
                for (int b = 0; b < n; b++) {
                    for (int i = 0; i < n; i++) pi[i] = (a * i + b) % n;
                    inverse(pi, sigma);
                    int I = toric(pi, sigma, &sum_rowmin, &tot, &bu, &bv);
                    checked++;
                    if (I > maxI) { maxI = I; memcpy(argI, pi, sizeof pi); }
                    if (sum_rowmin > maxRA) { maxRA = sum_rowmin; memcpy(argRA, pi, sizeof pi); }
                    if (I > bound(n)) { if (!nviol) memcpy(viol, pi, sizeof pi); nviol++; }
                    if (sum_rowmin > (long)n * bound(n)) nviolRA++;
                }
            }
            double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
            report("affine", checked, maxI, argI, maxRA, argRA, nviol, viol, secs);
            if (nviolRA) printf("{\"mode\":\"affine\",\"n\":%d,\"violations_RA\":%ld}\n", n, nviolRA);
        }
    } else if (!strcmp(mode, "rand") || !strcmp(mode, "climb")) {
        n = atoi(argv[2]);
        long trials = atol(argv[3]);
        rng_state = argc > 4 ? strtoull(argv[4], NULL, 10) : 12345;
        if (!rng_state) rng_state = 1;
        int maxI = -1; long maxRA = -1, nviol = 0, nviolRA = 0; long long checked = 0;
        memset(argI, 0, sizeof argI); memset(argRA, 0, sizeof argRA); memset(viol, 0, sizeof viol);
        int climb = !strcmp(mode, "climb");
        for (long t = 0; t < trials; t++) {
            for (int i = 0; i < n; i++) pi[i] = i;
            for (int i = n - 1; i > 0; i--) { int j = rnd() % (i + 1); int s = pi[i]; pi[i] = pi[j]; pi[j] = s; }
            inverse(pi, sigma);
            int I = toric(pi, sigma, &sum_rowmin, &tot, &bu, &bv);
            long RA = sum_rowmin;
            checked++;
            if (climb) {
                int improved = 1;
                while (improved) {
                    improved = 0;
                    for (int i = 0; i < n && !improved; i++)
                        for (int j = i + 1; j < n && !improved; j++) {
                            int s = pi[i]; pi[i] = pi[j]; pi[j] = s;
                            inverse(pi, sigma);
                            long sr, tt; int b1, b2;
                            int I2 = toric(pi, sigma, &sr, &tt, &b1, &b2);
                            checked++;
                            if (I2 > I || (I2 == I && sr > RA)) { I = I2; RA = sr; improved = 1; }
                            else { s = pi[i]; pi[i] = pi[j]; pi[j] = s; }
                        }
                }
                inverse(pi, sigma);
            }
            if (I > maxI) { maxI = I; memcpy(argI, pi, sizeof pi); }
            if (RA > maxRA) { maxRA = RA; memcpy(argRA, pi, sizeof pi); }
            if (I > bound(n)) { if (!nviol) memcpy(viol, pi, sizeof pi); nviol++; }
            if (RA > (long)n * bound(n)) nviolRA++;
        }
        double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
        report(mode, checked, maxI, argI, maxRA, argRA, nviol, viol, secs);
        printf("{\"mode\":\"%s\",\"n\":%d,\"violations_RA\":%ld}\n", mode, n, nviolRA);
    } else return 2;
    return 0;
}
