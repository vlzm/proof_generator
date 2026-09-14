/* toric_inversions.c — the purely combinatorial part of H13: I(pi), the minimum
 * number of inversions of the relabelled line over all n^2 double cuts
 * (equivalently over the toric class of pi: simultaneous cyclic rotation of
 * positions and of values).  Definitions follow experiments/line_model.py and
 * experiments/line_profile.c:
 *
 *     line(q, c):  w_j = pi[(q + 1 + j) mod n] - (q + 1 - c)   (mod n)
 *     I(pi)     =  min over the n^2 pairs (q, c) of inv(w).
 *
 * With alpha = q + 1 (position origin) and beta = q + 1 - c (value origin) the
 * pair (alpha, beta) runs over all of Z_n x Z_n, and
 *     w_j = pi[(alpha + j) mod n] - beta  (mod n),
 * which is the parametrisation used below.  The whole n x n table of inv is
 * built in O(n^2) from the two increment identities
 *     inv(alpha+1, beta) - inv(alpha, beta) = n - 1 - 2 * w_0,
 *     inv(alpha, beta+1) - inv(alpha, beta) = n - 1 - 2 * t_0,
 * where w_0 is the first value of the line and t_0 the line position of its
 * smallest value; mode "verify" checks the table against the direct definition.
 *
 * Modes (one JSON object on stdout per run):
 *   verify  n [samples]   table vs direct definition (all pi if n <= 8)
 *   exhaust n             all pi: max I, its argmax count, the averaging bound
 *                         I <= C(n,2) - T(pi)/n^2 and its slack, min of T
 *   free    n             all pi (n <= 10): I(pi) vs the minimum of the same
 *                         objective over all 2^(n-1) bipartitions (relaxation)
 *   climb   n restarts s  hill climbing on transpositions, maximise I
 *   affine  n             all pi(i) = a*i + b with gcd(a, n) = 1
 *
 * Version toric_inv-1.0.  Usage: gcc -O2 -o toric_inversions experiments/toric_inversions.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 40
static int n;
static int tbl[MAXN][MAXN];

static unsigned long long rs = 88172645463325252ULL;
static unsigned long long xr(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static void shuffle(int *a) { for (int i = n - 1; i > 0; i--) { int j = xr() % (i + 1); int t = a[i]; a[i] = a[j]; a[j] = t; } }

static int inv_direct(const int *pi, int alpha, int beta) {
    int w[MAXN], c = 0;
    for (int j = 0; j < n; j++) { int v = (pi[(alpha + j) % n] - beta) % n; if (v < 0) v += n; w[j] = v; }
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) if (w[i] > w[j]) c++;
    return c;
}

/* fills tbl[alpha][beta] and returns I(pi) */
static int compute_I(const int *pi) {
    int ipi[MAXN], col0[MAXN];
    for (int i = 0; i < n; i++) ipi[pi[i]] = i;
    col0[0] = inv_direct(pi, 0, 0);
    for (int a = 0; a + 1 < n; a++) col0[a + 1] = col0[a] + n - 1 - 2 * pi[a];
    int best = 1 << 30;
    for (int a = 0; a < n; a++) {
        int cur = col0[a];
        for (int b = 0; b < n; b++) {
            tbl[a][b] = cur;
            if (cur < best) best = cur;
            int t = (ipi[b] - a) % n; if (t < 0) t += n;
            cur += n - 1 - 2 * t;
        }
    }
    return best;
}

/* T(pi) = sum over ordered pairs a != b of ((x_b - x_a) mod n) * ((y_b - y_a) mod n)
 * = sum over all n^2 cuts of the number of concordant (non-inverted) pairs. */
static long long Tvalue(const int *pi) {
    long long T = 0;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j) {
        int p = (j - i) % n; if (p < 0) p += n;
        int v = (pi[j] - pi[i]) % n; if (v < 0) v += n;
        T += (long long)p * v;
    }
    return T;
}

static int next_perm(int *a, int m) {
    int i = m - 2; while (i >= 0 && a[i] >= a[i + 1]) i--;
    if (i < 0) return 0;
    int j = m - 1; while (a[j] <= a[i]) j--;
    int t = a[i]; a[i] = a[j]; a[j] = t;
    for (int l = i + 1, r = m - 1; l < r; l++, r--) { t = a[l]; a[l] = a[r]; a[r] = t; }
    return 1;
}

/* minimum over all 2^n sign vectors x of  (C(n,2) + sum_{a<b} eps_ab x_a x_b) / 2,
 * eps_ab = +1 if the pair is inverted at the base cut (0,0), -1 otherwise.
 * This is the relaxation of I(pi): quadrant colourings are replaced by all
 * bipartitions (|E symmetric-difference cut(U)| over all U). */
static int free_min(const int *pi) {
    int eps[MAXN][MAXN];
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) eps[i][j] = (i < j) ? (pi[i] > pi[j] ? 1 : -1) : 0;
    int best = 1 << 30;
    for (long U = 0; U < (1L << (n - 1)); U++) {
        int x[MAXN];
        for (int i = 0; i < n; i++) x[i] = ((U >> i) & 1) ? -1 : 1;
        int R = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) R += eps[i][j] * x[i] * x[j];
        int val = (n * (n - 1) / 2 + R) / 2;
        if (val < best) best = val;
    }
    return best;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s mode n [args]\n", argv[0]); return 2; }
    const char *mode = argv[1];
    n = atoi(argv[2]);
    if (n < 2 || n > MAXN) return 2;
    int tgt = ((n - 1) * (n - 1)) / 4;
    long long C2 = (long long)n * (n - 1) / 2;

    if (!strcmp(mode, "verify")) {
        long samples = argc > 3 ? atol(argv[3]) : 200;
        long bad = 0, checked = 0;
        int pi[MAXN];
        if (n <= 8) {
            for (int i = 0; i < n; i++) pi[i] = i;
            do {
                compute_I(pi);
                for (int a = 0; a < n; a++) for (int b = 0; b < n; b++)
                    if (tbl[a][b] != inv_direct(pi, a, b)) bad++;
                checked++;
            } while (next_perm(pi, n));
        } else {
            for (int i = 0; i < n; i++) pi[i] = i;
            for (long s = 0; s < samples; s++) {
                shuffle(pi); compute_I(pi);
                for (int a = 0; a < n; a++) for (int b = 0; b < n; b++)
                    if (tbl[a][b] != inv_direct(pi, a, b)) bad++;
                checked++;
            }
        }
        printf("{\"mode\": \"verify\", \"n\": %d, \"perms_checked\": %ld, \"mismatches\": %ld}\n", n, checked, bad);
        return bad ? 1 : 0;
    }

    if (!strcmp(mode, "exhaust")) {
        int pi[MAXN]; for (int i = 0; i < n; i++) pi[i] = i;
        int maxI = -1; long argcnt = 0; long total = 0;
        long long minT = -1; int bound_violations = 0; long long worst_slack = -1;
        int maxI_pi[MAXN];
        do {
            int I = compute_I(pi);
            long long T = Tvalue(pi);
            total++;
            if (minT < 0 || T < minT) minT = T;
            /* averaging bound: I <= C(n,2) - T/n^2, i.e. n^2 * (C2 - I) >= T */
            if ((long long)n * n * (C2 - I) < T) bound_violations++;
            long long slack = (C2 - (long long)I) * n * n - T; /* >= 0 */
            if (worst_slack < 0 || slack < worst_slack) worst_slack = slack;
            if (I > maxI) { maxI = I; argcnt = 1; memcpy(maxI_pi, pi, sizeof(int) * n); }
            else if (I == maxI) argcnt++;
        } while (next_perm(pi, n));
        printf("{\"mode\": \"exhaust\", \"n\": %d, \"perms\": %ld, \"target_floor_(n-1)^2/4\": %d, "
               "\"max_I\": %d, \"argmax_count\": %ld, \"first_argmax\": [", n, total, tgt, maxI, argcnt);
        for (int i = 0; i < n; i++) printf("%s%d", i ? ", " : "", maxI_pi[i]);
        printf("], \"min_T\": %lld, \"rearrangement_bound_n^2(n^2-1)/6\": %lld, "
               "\"averaging_bound_violations\": %d, \"min_averaging_slack\": %lld, "
               "\"averaging_bound_value_floor\": %lld}\n",
               minT, (long long)n * n * ((long long)n * n - 1) / 6, bound_violations, worst_slack,
               (long long)((C2 * n * n - (long long)n * n * ((long long)n * n - 1) / 6) / ((long long)n * n)));
        return 0;
    }

    if (!strcmp(mode, "free")) {
        int pi[MAXN]; for (int i = 0; i < n; i++) pi[i] = i;
        int maxfree = -1, maxgap = 0, maxI = -1; long ndiff = 0, total = 0;
        int gap_pi[MAXN];
        do {
            int I = compute_I(pi), F = free_min(pi);
            total++;
            if (F > maxfree) maxfree = F;
            if (I > maxI) maxI = I;
            if (I != F) ndiff++;
            if (I - F > maxgap) { maxgap = I - F; memcpy(gap_pi, pi, sizeof(int) * n); }
        } while (next_perm(pi, n));
        printf("{\"mode\": \"free\", \"n\": %d, \"perms\": %ld, \"target\": %d, \"max_I\": %d, "
               "\"max_free_min\": %d, \"perms_with_I_gt_free\": %ld, \"max_gap\": %d, \"gap_pi\": [",
               n, total, tgt, maxI, maxfree, ndiff, maxgap);
        for (int i = 0; i < n; i++) printf("%s%d", i ? ", " : "", gap_pi[i]);
        printf("]}\n");
        return 0;
    }

    if (!strcmp(mode, "climb")) {
        long restarts = argc > 3 ? atol(argv[3]) : 100;
        if (argc > 4) rs = (unsigned long long)atoll(argv[4]) * 2862933555777941757ULL + 3037000493ULL;
        int pi[MAXN], bp[MAXN]; int gbest = -1;
        for (long it = 0; it < restarts; it++) {
            for (int i = 0; i < n; i++) pi[i] = i;
            shuffle(pi);
            int cur = compute_I(pi), improved = 1;
            while (improved) {
                improved = 0;
                for (int i = 0; i < n && !improved; i++) for (int j = i + 1; j < n; j++) {
                    int t = pi[i]; pi[i] = pi[j]; pi[j] = t;
                    int v = compute_I(pi);
                    if (v > cur) { cur = v; improved = 1; break; }
                    t = pi[i]; pi[i] = pi[j]; pi[j] = t;
                }
            }
            if (cur > gbest) { gbest = cur; memcpy(bp, pi, sizeof(int) * n); }
        }
        printf("{\"mode\": \"climb\", \"n\": %d, \"restarts\": %ld, \"target\": %d, \"max_I\": %d, "
               "\"excess\": %d, \"argmax\": [", n, restarts, tgt, gbest, gbest - tgt);
        for (int i = 0; i < n; i++) printf("%s%d", i ? ", " : "", bp[i]);
        printf("]}\n");
        return 0;
    }

    if (!strcmp(mode, "refl")) {
        /* the neighbourhood of the maximisers: every reflection pi_h(i) = h - i mod n
         * with one or two transpositions of values applied */
        int pi[MAXN], base[MAXN]; int best1 = -1, best2 = -1; int b1[MAXN], b2[MAXN];
        long cnt = 0;
        for (int h = 0; h < n; h++) {
            for (int i = 0; i < n; i++) base[i] = ((h - i) % n + n) % n;
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
                memcpy(pi, base, sizeof(int) * n);
                int t = pi[i]; pi[i] = pi[j]; pi[j] = t;
                int I = compute_I(pi); cnt++;
                if (I > best1) { best1 = I; memcpy(b1, pi, sizeof(int) * n); }
                for (int k = 0; k < n; k++) for (int l = k + 1; l < n; l++) {
                    int u = pi[k]; pi[k] = pi[l]; pi[l] = u;
                    int J = compute_I(pi); cnt++;
                    if (J > best2) { best2 = J; memcpy(b2, pi, sizeof(int) * n); }
                    u = pi[k]; pi[k] = pi[l]; pi[l] = u;
                }
            }
        }
        printf("{\"mode\": \"refl\", \"n\": %d, \"target\": %d, \"tested\": %ld, "
               "\"max_I_one_transposition\": %d, \"max_I_two_transpositions\": %d, \"excess\": %d, \"argmax2\": [",
               n, tgt, cnt, best1, best2, (best1 > best2 ? best1 : best2) - tgt);
        for (int i = 0; i < n; i++) printf("%s%d", i ? ", " : "", b2[i]);
        printf("]}\n");
        return 0;
    }

    if (!strcmp(mode, "affine")) {
        int pi[MAXN]; int best = -1, ba = 0, bb = 0;
        for (int a = 1; a < n; a++) {
            int x = a, y = n; while (y) { int t = x % y; x = y; y = t; }
            if (x != 1) continue;
            for (int b = 0; b < n; b++) {
                for (int i = 0; i < n; i++) pi[i] = (a * i + b) % n;
                int I = compute_I(pi);
                if (I > best) { best = I; ba = a; bb = b; }
            }
        }
        printf("{\"mode\": \"affine\", \"n\": %d, \"target\": %d, \"max_I\": %d, \"excess\": %d, "
               "\"a\": %d, \"b\": %d}\n", n, tgt, best, best - tgt, ba, bb);
        return 0;
    }
    fprintf(stderr, "unknown mode %s\n", mode);
    return 2;
}
