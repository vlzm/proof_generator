/* search_C27.c -- adversarial search for permutations that make the
 * carrier-route expression of construction N1 large:
 *
 *     maximize over pi in S_n   min_c [ 2 F_c - S_c + R_c ] - B_n .
 *
 * H10 / C27 says this is always <= 1.  Exhaustive enumeration is only feasible
 * up to n = 11 (checks/check_C27_fast.c); for larger n this does randomized
 * hill climbing over transposition neighbourhoods, which is a search for
 * counterexamples, not a proof (AGENTS.md rule 9).
 *
 * Build: gcc -O2 -o experiments/search_C27 experiments/search_C27.c
 * Run:   ./experiments/search_C27 <n> <restarts> <seed>
 * Prints the best value found, how often the maximum was reached and one
 * witness per distinct best value; SAMPLED coverage.
 */

#include <string.h>
#include <time.h>
#include "../checks/carrier_core.h"

static unsigned long long rng_state;
static unsigned long long rnd(void) {         /* xorshift64* */
    unsigned long long x = rng_state;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    rng_state = x;
    return x * 2685821657736338717ULL;
}
static int rndint(int m) { return (int)(rnd() % (unsigned)m); }

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s n restarts seed\n", argv[0]); return 1; }
    int n = atoi(argv[1]);
    long restarts = atol(argv[2]);
    rng_state = (unsigned long long)atoll(argv[3]) * 2862933555777941757ULL + 3037000493ULL;
    if (n < 3 || n > MAXN) { fprintf(stderr, "n out of range\n"); return 1; }
    cc_set_n(n);
    int Bn = n * (n - 1) / 2;
    int pi[MAXN], bestpi[MAXN];
    int global = -100000, hits = 0;
    clock_t t0 = clock();

    /* structured family first: all affine permutations pi(i) = a*i + b,
       gcd(a,n) = 1 (a = n-1 are the reflections) */
    int affbest = -100000, affa = 0, affb = 0;
    for (int a = 1; a < n; a++) {
        int g = a, h = n;
        while (h) { int t = g % h; g = h; h = t; }
        if (g != 1) continue;
        for (int b = 0; b < n; b++) {
            for (int i = 0; i < n; i++) pi[i] = (a * i + b) % n;
            int v = cc_min_over_c(pi, NULL, NULL);
            if (v > affbest) { affbest = v; affa = a; affb = b; }
        }
    }
    printf("affine family: max min_c = %d (= B_n %+d) at pi(i) = %d*i + %d\n",
           affbest, affbest - Bn, affa, affb);
    for (int i = 0; i < n; i++) bestpi[i] = (affa * i + affb) % n;
    global = affbest; hits = 1;

    for (long r = 0; r < restarts; r++) {
        for (int i = 0; i < n; i++) pi[i] = i;
        for (int i = n - 1; i > 0; i--) { int j = rndint(i + 1); int t = pi[i]; pi[i] = pi[j]; pi[j] = t; }
        int cur = cc_min_over_c(pi, NULL, NULL);
        int improved = 1;
        while (improved) {                       /* steepest ascent on swaps */
            improved = 0;
            int bi = -1, bj = -1, bv = cur;
            for (int i = 0; i < n; i++)
                for (int j = i + 1; j < n; j++) {
                    int t = pi[i]; pi[i] = pi[j]; pi[j] = t;
                    int v = cc_min_over_c(pi, NULL, NULL);
                    t = pi[i]; pi[i] = pi[j]; pi[j] = t;
                    if (v > bv) { bv = v; bi = i; bj = j; }
                }
            if (bi >= 0) { int t = pi[bi]; pi[bi] = pi[bj]; pi[bj] = t; cur = bv; improved = 1; }
        }
        if (cur > global) {
            global = cur; hits = 1; memcpy(bestpi, pi, sizeof(int) * n);
            printf("new best %d (= B_n %+d) at pi =", cur, cur - Bn);
            for (int i = 0; i < n; i++) printf(" %d", pi[i]);
            printf("\n"); fflush(stdout);
        } else if (cur == global) hits++;
    }
    printf("n=%d B_n=%d restarts=%ld best=%d excess=%+d hits=%d time=%.1fs SAMPLED\n",
           n, Bn, restarts, global, global - Bn, hits, (double)(clock() - t0) / CLOCKS_PER_SEC);
    printf("witness:");
    for (int i = 0; i < n; i++) printf(" %d", bestpi[i]);
    printf("\n");
    return 0;
}
