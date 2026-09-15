/* h13i_avg.c -- exhaustive check of a candidate route to H13-I (PLAN Sec 8):
 *
 *   I(pi) = min_{p0,v0} inv <= floor((n-1)^2/4)   for every permutation pi of Z_n.
 *
 * Candidate lemma (new this session, not previously tried per
 * docs/notes/h13_line_model.md Sec 6): average the ADAPTIVE minimum over one
 * coordinate only, instead of the uniform average over all n^2 cuts (which is
 * already known to fail, PLAN Sec 8 / h13_line_model.md Sec 5).
 *
 * For a fixed position cut p0 let
 *   g(p0) = min_{v0 in Z_n} inv( (pi(p0+j) - v0) mod n : j = 0..n-1 )
 * (the best value cut for that position cut).  Trivially I(pi) <= g(p0) for
 * every p0, so it suffices that SOME p0 has g(p0) <= floor((n-1)^2/4); by
 * pigeonhole this holds if the AVERAGE of g(p0) over p0 = 0..n-1 does:
 *
 *   CANDIDATE LEMMA:  sum_{p0=0}^{n-1} g(p0)  <=  n * floor((n-1)^2/4).
 *
 * g(p0) is computed via the exact telescoping recurrence
 *   inv(p0, v+1) - inv(p0, v) = n - 1 - 2*pos(v),   pos(v) = position (in the
 * p0-rotated line) of value v -- i.e. moving the value v from being the
 * smallest label to being the largest, all other relative order unchanged.
 * (Derivation and a hand check are in docs/notes/h13_line_model.md Sec 7.)
 * This turns each g(p0) into an O(n) walk after one O(n^2) inversion count,
 * so the whole check is O(n^3) per permutation -- exhaustive up to n = 12 is
 * feasible (n = 12: 479001600 permutations).
 *
 * This program only reports numbers; it does NOT constitute a proof (rule 5,
 * AGENTS.md). A finite check never proves a general claim.
 *
 * Usage: h13i_avg n
 * Version h13i_avg-1.0 (session 9).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXN 13

static int inv_count(const int *seq, int n) {
    int inv = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (seq[i] > seq[j]) inv++;
    return inv;
}

static int g_of_p0(const int *line, int n, int *pos_buf) {
    int inv0 = inv_count(line, n);
    for (int idx = 0; idx < n; idx++) pos_buf[line[idx]] = idx;
    int S = 0, best = inv0;
    for (int v = 0; v < n - 1; v++) {
        int D = n - 1 - 2 * pos_buf[v];
        S += D;
        int cur = inv0 + S;
        if (cur < best) best = cur;
    }
    return best;
}

static int perm[MAXN], line[MAXN], posbuf[MAXN], worst_perm[MAXN];

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s n\n", argv[0]); return 2; }
    int n = atoi(argv[1]);
    if (n < 2 || n > MAXN) { fprintf(stderr, "n out of range 2..%d\n", MAXN); return 2; }
    for (int i = 0; i < n; i++) perm[i] = i;
    long worst_sum = -1, count = 0;
    int bnd = (n - 1) * (n - 1) / 4;
    clock_t t0 = clock();
    do {
        long tot = 0;
        for (int p0 = 0; p0 < n; p0++) {
            for (int j = 0; j < n; j++) line[j] = perm[(p0 + j) % n];
            tot += g_of_p0(line, n, posbuf);
        }
        if (tot > worst_sum) { worst_sum = tot; memcpy(worst_perm, perm, sizeof(int) * n); }
        count++;
        int i = n - 2;
        while (i >= 0 && perm[i] >= perm[i + 1]) i--;
        if (i < 0) break;
        int j = n - 1;
        while (perm[j] <= perm[i]) j--;
        int t = perm[i]; perm[i] = perm[j]; perm[j] = t;
        for (int l = i + 1, r = n - 1; l < r; l++, r--) { t = perm[l]; perm[l] = perm[r]; perm[r] = t; }
    } while (1);
    double dt = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("{\"version\":\"h13i_avg-1.0\",\"n\":%d,\"bound_per_p0\":%d,\"bound_times_n\":%d,"
           "\"worst_sum\":%ld,\"worst_avg\":%.6f,\"count\":%ld,\"time_s\":%.2f,\"status\":\"%s\","
           "\"worst_pi\":[",
           n, bnd, bnd * n, worst_sum, (double)worst_sum / n, count, dt,
           worst_sum <= (long)bnd * n ? "OK" : "VIOLATION");
    for (int i = 0; i < n; i++) printf("%d%s", worst_perm[i], i + 1 < n ? "," : "");
    printf("]}\n");
    return 0;
}
