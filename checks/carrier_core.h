/* carrier_core.h -- the quantity  2F_c - S_c + R_c  of the carrier-route variant
 * of construction N1, reimplemented from the specification (PROBLEM.md 5.2,
 * N1 sections 1-3) independently of constructions/strict_upper.py.
 *
 * Conventions (AGENTS.md rule 16): signed shortest step d with -n/2 < d <= n/2,
 * so an antipodal step is positive.  Per nontrivial cycle C of f_c(i)=pi(i)+c:
 *   F(C) = sum |d_j|,  M(C) = max edge load,  t(C) = #{j : d_{j-1}<0<d_j},
 *   E(C) = |C| - 2t(C),  S(C) = 2M(C) + E(C) - 1,
 * carrier b(C) = first vertex of the cycle (in f-order from its smallest
 * element) when all steps have one sign, else the smallest v with load(v) = M
 * and load(v-1) < M.  R_c = length of a shortest walk on Z_n from 0 to c that
 * visits every carrier.
 *
 * Cross-checked against strict_upper-1.0 by checks/check_C27.py.
 */

#ifndef CARRIER_CORE_H
#define CARRIER_CORE_H

#include <stdio.h>
#include <stdlib.h>

#define MAXN 60        /* carrier sets are kept as 64-bit bitmasks */
typedef unsigned long long cc_mask;

static int CC_n;

static inline void cc_set_n(int n) { CC_n = n; }

static inline int cc_signed_step(int a, int b) {
    int n = CC_n;
    int d = (b - a) % n;
    if (d < 0) d += n;
    if (2 * d > n) d -= n;
    return d;
}

/* shortest walk on Z_n from 0 to c visiting every vertex of reqmask.
 * The walk covers a contiguous arc; enumerate the gap left uncovered. */
static int cc_walk_cost(cc_mask reqmask, int c) {
    int n = CC_n;
    int pts[MAXN + 2], m = 0, v;
    reqmask |= ((cc_mask)1 << 0) | ((cc_mask)1 << c);
    for (v = 0; v < n; v++) if (reqmask & ((cc_mask)1 << v)) pts[m++] = v;
    int best = -1;
    for (int g = 0; g < m; g++) {
        int lo = pts[(g + 1) % m];
        int hi = pts[g];
        int L = (m == 1) ? 0 : (((hi - lo) % n) + n) % n;
        int x0 = (((0 - lo) % n) + n) % n;
        int xc = (((c - lo) % n) + n) % n;
        if (x0 > L || xc > L) { fprintf(stderr, "cc_walk_cost: arc error\n"); exit(2); }
        int costA = x0 + L + (L - xc);          /* down to lo, up to hi, back to c */
        int costB = (L - x0) + L + xc;          /* up to hi, down to lo, up to c   */
        if (best < 0 || costA < best) best = costA;
        if (costB < best) best = costB;
    }
    int up = n + c;                              /* full turn +, then on to c */
    int dn = n + (n - c) % n;                    /* full turn -, then on to c */
    if (up < best) best = up;
    if (dn < best) best = dn;
    return best;
}

/* 2F_c - S_c + R_c for one shift; optional outputs F, S, R, carrier mask. */
static int cc_shift_value(const int *pi, int c, int *out_F, int *out_S, int *out_R, cc_mask *out_mask) {
    int n = CC_n;
    int f[MAXN], seen[MAXN], cyc[MAXN], steps[MAXN], loads[MAXN];
    int i, j;
    for (i = 0; i < n; i++) { f[i] = (pi[i] + c) % n; seen[i] = 0; }
    int F = 0, S = 0;
    cc_mask mask = 0;
    for (i = 0; i < n; i++) {
        if (seen[i] || f[i] == i) { seen[i] = 1; continue; }
        int k = 0, x = i;
        while (!seen[x]) { seen[x] = 1; cyc[k++] = x; x = f[x]; }
        for (j = 0; j < n; j++) loads[j] = 0;
        int Fc = 0, allpos = 1, allneg = 1;
        for (j = 0; j < k; j++) {
            int d = cc_signed_step(cyc[j], cyc[(j + 1) % k]);
            steps[j] = d;
            Fc += (d > 0 ? d : -d);
            if (d > 0) { allneg = 0; for (int q = 0; q < d; q++) loads[(cyc[j] + q) % n]++; }
            else       { allpos = 0; for (int q = 0; q < -d; q++) loads[(((cyc[j] - 1 - q) % n) + n) % n]++; }
        }
        int M = 0;
        for (j = 0; j < n; j++) if (loads[j] > M) M = loads[j];
        int t = 0;
        for (j = 0; j < k; j++) if (steps[(j + k - 1) % k] < 0 && steps[j] > 0) t++;
        int Sc = 2 * M + (k - 2 * t) - 1;
        int b;
        if (allpos || allneg) {
            b = cyc[0];
        } else {
            b = -1;
            for (j = 0; j < n; j++) if (loads[j] == M && loads[(j + n - 1) % n] < M) { b = j; break; }
            if (b < 0) { fprintf(stderr, "cc: carrier not found\n"); exit(2); }
        }
        if (mask & ((cc_mask)1 << b)) { fprintf(stderr, "cc: duplicate carrier\n"); exit(2); }
        mask |= (cc_mask)1 << b;
        F += Fc; S += Sc;
    }
    int R = cc_walk_cost(mask, c);
    if (out_F) *out_F = F;
    if (out_S) *out_S = S;
    if (out_R) *out_R = R;
    if (out_mask) *out_mask = mask;
    return 2 * F - S + R;
}

/* min over all shifts c; also returns the sum over c if wanted. */
static int cc_min_over_c(const int *pi, long long *out_sum, int *out_argmin) {
    int best = -1, arg = -1;
    long long sum = 0;
    for (int c = 0; c < CC_n; c++) {
        int v = cc_shift_value(pi, c, NULL, NULL, NULL, NULL);
        sum += v;
        if (best < 0 || v < best) { best = v; arg = c; }
    }
    if (out_sum) *out_sum = sum;
    if (out_argmin) *out_argmin = arg;
    return best;
}

#endif
