/* h13i_search.c — adversarial search for a counterexample to H13-I:
 * (h13i_search-1.1: simulated annealing + reflection-seeded restarts,
 * see docs/notes/h13i_verdict.md for why plain hill-climbing (1.0) was too
 * weak a probe to trust a negative result)
 * claim I(pi) <= floor((n-1)^2/4) for every permutation pi of Z_n, where
 * I(pi) = min over double cuts (q, c) of inv(w) (see docs/notes/h13_line_model.md
 * §0; equivalent definition here: min over position-rotation a and
 * value-shift b of inversions of the sequence row[j] = (pi[(a+j)%n]-b)%n).
 *
 * Exhaustive check already covers 4<=n<=10 (experiments/line_profile.c,
 * session 8). This tool does NOT try to be exhaustive: it computes I(pi)
 * for a single permutation in O(n^2) (fast incremental recurrences, no BFS,
 * no table lookup) and runs a randomized hill-climbing / restart search that
 * tries to MAXIMIZE I(pi) - floor((n-1)^2/4) for n beyond the exhaustive
 * range, to probe whether the bound can be violated. A negative result here
 * is evidence, not a proof, and is reported as such (AGENTS.md rule 9).
 *
 * O(n^2) I(pi) computation:
 *  - pos-shift recurrence (fixed b=0): inv(a+1,0) = inv(a,0) + (n-1-2*pi[a]).
 *  - value-shift recurrence (fixed a): inv(a,b+1) = inv(a,b) + (n-1-2*pos_a[b]),
 *    where pos_a[v] = (piinv[v]-a) mod n is the row-position of value v in
 *    the rotation starting at a.
 * inv(0,0) computed once with an O(n^2) double loop (n is small enough that
 * this dominates nothing at the search sizes used here).
 *
 * Usage:
 *   h13i_search n seed iters [restarts]
 * Prints one line per new best found: n, iters-so-far, I, bound, margin(I-bound), pi.
 * Version h13i_search-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int n;
static int pi_arr[4096], piinv[4096], pos_a[4096];

static long floor_bound(int n) {
    long m = n - 1;
    return (m * m) / 4;
}

static long inv0_full(const int *row) {
    long c = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (row[i] > row[j]) c++;
    return c;
}

/* returns I(pi) = min over (a,b) */
static long I_of_pi(void) {
    for (int v = 0; v < n; v++) piinv[pi_arr[v]] = v;
    int row0[4096];
    for (int j = 0; j < n; j++) row0[j] = pi_arr[j];
    long inv_a0 = inv0_full(row0); /* inv(a=0, b=0) */
    long best = inv_a0;
    long inv_a = inv_a0;
    for (int a = 0; a < n; a++) {
        for (int v = 0; v < n; v++) {
            int p = piinv[v] - a; if (p < 0) p += n;
            pos_a[v] = p;
        }
        long inv_b = inv_a;
        if (inv_b < best) best = inv_b;
        for (int b = 0; b < n - 1; b++) {
            inv_b += (long)(n - 1) - 2 * pos_a[b];
            if (inv_b < best) best = inv_b;
        }
        /* advance a: inv(a+1,0) = inv(a,0) + (n-1-2*pi_arr[a]) */
        inv_a += (long)(n - 1) - 2 * pi_arr[a];
    }
    return best;
}

static unsigned long rng_state;
static unsigned long xrand(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}

static double xrand01(void) { return (double)(xrand() % 1000000007UL) / 1000000007.0; }

/* seed[0]: 0 = uniform random permutation, 1 = reflection pi(i) = (h-i) mod n
 * for random h, 2 = identity (edge case seed) */
static void seed_pi(int mode) {
    if (mode == 1) {
        int h = xrand() % n;
        for (int i = 0; i < n; i++) { int v = h - i; v %= n; if (v < 0) v += n; pi_arr[i] = v; }
    } else if (mode == 2) {
        for (int i = 0; i < n; i++) pi_arr[i] = i;
    } else {
        for (int i = 0; i < n; i++) pi_arr[i] = i;
        for (int i = n - 1; i > 0; i--) {
            int j = xrand() % (i + 1);
            int t = pi_arr[i]; pi_arr[i] = pi_arr[j]; pi_arr[j] = t;
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s n seed iters [restarts]\n", argv[0]); return 2; }
    n = atoi(argv[1]);
    unsigned long seed = strtoul(argv[2], NULL, 10);
    long iters = atol(argv[3]);
    int restarts = argc >= 5 ? atoi(argv[4]) : 1;
    if (n < 4 || n > 4000) { fprintf(stderr, "n out of range\n"); return 2; }
    rng_state = seed ? seed : 88172645463325252UL;

    long bound = floor_bound(n);
    long global_best = -1;
    int best_pi[4096];

    for (int r = 0; r < restarts; r++) {
        int mode = r % 3; /* cycle: random, reflection, random, reflection, ... */
        seed_pi(mode);
        long cur = I_of_pi();
        double T0 = (double)n / 2.0, T;
        for (long it = 0; it < iters; it++) {
            T = T0 * (1.0 - (double)it / (double)iters) + 0.01;
            int i = xrand() % n, j = xrand() % n;
            if (i == j) continue;
            int t = pi_arr[i]; pi_arr[i] = pi_arr[j]; pi_arr[j] = t;
            long cand = I_of_pi();
            long delta = cand - cur;
            int accept = (delta >= 0) || (xrand01() < __builtin_exp((double)delta / T));
            if (accept) {
                cur = cand;
            } else {
                t = pi_arr[i]; pi_arr[i] = pi_arr[j]; pi_arr[j] = t;
            }
            if (cur > global_best) {
                global_best = cur;
                memcpy(best_pi, pi_arr, sizeof(int) * n);
                printf("n=%d restart=%d mode=%d iter=%ld I=%ld bound=%ld margin=%ld pi=[", n, r, mode, it, cur, bound, cur - bound);
                for (int k = 0; k < n; k++) printf("%d%s", best_pi[k], k + 1 < n ? "," : "");
                printf("]\n");
                fflush(stdout);
            }
        }
    }
    printf("FINAL n=%d bound=%ld best_I=%ld margin=%ld\n", n, bound, global_best, global_best - bound);
    return 0;
}
