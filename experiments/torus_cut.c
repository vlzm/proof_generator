/* torus_cut.c — fast evaluation of I(pi) = min_{q,c} inv(w_{q,c}) (H13-I,
 * docs/notes/h13_line_model.md), using the closed-form / increment lemmas of
 * docs/proofs/C37_cut_inversion_formula.md instead of the O(n^2) per-cut
 * brute force of experiments/line_profile.c. This drops the cost of scanning
 * all n^2 cuts of one permutation from O(n^4) to O(n^2), which is what makes
 * the exhaustive checks below reach n = 11..13 and the structured/random
 * checks reach n up to a few thousand within a session budget.
 *
 * Parametrization: a = q + 1 (window start), b = q + 1 - c (value shift);
 * w_j = (pi[(a+j) mod n] - b) mod n, j = 0..n-1; f(a,b) = inv(w). Equivalent
 * to line_profile.c's (q, c) via a = q+1 mod n, c = a - b mod n.
 *
 * Modes:
 *   exhaustive N            all N! permutations (Heap's algorithm); tracks
 *                           max_pi I(pi) and reports whether it exceeds
 *                           floor((N-1)^2/4).
 *   family N [step]         reflections, the coprime-slope affine family,
 *                           and random samples at size N (see set_reflection,
 *                           set_affine, set_random below); step subsamples
 *                           b for the affine family when N is large.
 *   search N iters restarts seed
 *                           simulated annealing over S_N (transpositions and
 *                           segment reversals) trying to MAXIMIZE I(pi), i.e.
 *                           trying to refute I(pi) <= floor((n-1)^2/4).
 *
 * Usage: gcc -O3 -o torus_cut torus_cut.c && ./torus_cut exhaustive 11
 * Version torus_cut-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

static int n;
static int *pi_, *inv_pi, *a0_vals;

/* f(a, 0) for a = 0..n-1 via the increment lemma f(a+1,0) - f(a,0) =
 * n - 1 - 2*pi(a), then f(a, b) for b = 0..n-1 via f(a,b+1) - f(a,b) =
 * n - 1 - 2*((pi^{-1}(b) - a) mod n). Returns min over all n^2 cuts. */
static long compute_I(void) {
    for (int i = 0; i < n; i++) inv_pi[pi_[i]] = i;
    long f00 = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (pi_[i] > pi_[j]) f00++;
    a0_vals[0] = f00;
    long cur = f00;
    for (int a = 0; a < n - 1; a++) {
        cur += (n - 1 - 2 * pi_[a]);
        a0_vals[a + 1] = cur;
    }
    long best = f00;
    for (int a = 0; a < n; a++) {
        long g = a0_vals[a];
        if (g < best) best = g;
        for (int b = 0; b < n - 1; b++) {
            g += (n - 1 - 2 * ((inv_pi[b] - a + n) % n));
            if (g < best) best = g;
        }
    }
    return best;
}

/* ---- exhaustive mode ---- */

static long total_perms = 0;
static long best_I_seen = -1;

static void process(void) {
    total_perms++;
    long I = compute_I();
    if (I > best_I_seen) {
        best_I_seen = I;
        printf("new best I=%ld at perm:", I);
        for (int i = 0; i < n; i++) printf(" %d", pi_[i]);
        printf("\n");
        fflush(stdout);
    }
}

static void heaps(int k) {
    if (k == 1) { process(); return; }
    for (int i = 0; i < k; i++) {
        heaps(k - 1);
        if (k % 2 == 0) { int t = pi_[i]; pi_[i] = pi_[k - 1]; pi_[k - 1] = t; }
        else { int t = pi_[0]; pi_[0] = pi_[k - 1]; pi_[k - 1] = t; }
    }
}

static int run_exhaustive(void) {
    for (int i = 0; i < n; i++) pi_[i] = i;
    heaps(n);
    long bound = (long)(n - 1) * (n - 1) / 4;
    printf("n=%d total_perms=%ld bound=%ld best_I=%ld %s\n", n, total_perms, bound, best_I_seen,
           best_I_seen <= bound ? "OK" : "VIOLATION");
    return best_I_seen <= bound ? 0 : 1;
}

/* ---- family mode ---- */

static long egcd(long a, long b, long *x, long *y) {
    if (b == 0) { *x = 1; *y = 0; return a; }
    long x1, y1;
    long g = egcd(b, a % b, &x1, &y1);
    *x = y1; *y = x1 - (a / b) * y1;
    return g;
}

static void set_reflection(int h) { for (int i = 0; i < n; i++) pi_[i] = ((h - i) % n + n) % n; }

static void set_affine(long amul, long b) {
    for (int i = 0; i < n; i++) pi_[i] = (int)(((amul * i + b) % n + n) % n);
}

static void set_random(unsigned seed) {
    for (int i = 0; i < n; i++) pi_[i] = i;
    srand(seed);
    for (int i = n - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        int t = pi_[i]; pi_[i] = pi_[j]; pi_[j] = t;
    }
}

static int run_family(int step_arg) {
    long bound = (long)(n - 1) * (n - 1) / 4;
    int ok = 1;

    long max_refl = -1;
    int hstep = step_arg > 0 ? step_arg : (n > 2000 ? n / 500 + 1 : 1);
    for (int h = 0; h < n; h += hstep) { set_reflection(h); long I = compute_I(); if (I > max_refl) max_refl = I; }
    printf("n=%d bound=%ld reflections(step=%d) max_I=%ld %s\n", n, bound, hstep, max_refl,
           max_refl <= bound ? "OK" : "VIOLATION");
    ok &= (max_refl <= bound);

    long max_affine = -1; long best_a = -1, best_b = -1;
    long bstep = n > 200 ? n / 50 + 1 : 1;
    for (long am = 1; am < n; am++) {
        long x, y; if (egcd(am, n, &x, &y) != 1) continue;
        for (long bb = 0; bb < n; bb += bstep) {
            set_affine(am, bb);
            long I = compute_I();
            if (I > max_affine) { max_affine = I; best_a = am; best_b = bb; }
        }
    }
    printf("n=%d bound=%ld affine max_I=%ld at a=%ld b=%ld %s\n", n, bound, max_affine, best_a, best_b,
           max_affine <= bound ? "OK" : "VIOLATION");
    ok &= (max_affine <= bound);

    long max_rand = -1;
    int nsamples = n > 500 ? 200 : 2000;
    for (int s = 0; s < nsamples; s++) { set_random(1000 + s); long I = compute_I(); if (I > max_rand) max_rand = I; }
    printf("n=%d bound=%ld random(%d samples) max_I=%ld %s\n", n, bound, nsamples, max_rand,
           max_rand <= bound ? "OK" : "VIOLATION");
    ok &= (max_rand <= bound);

    return ok ? 0 : 1;
}

/* ---- simulated-annealing search for a counterexample (maximize I) ---- */

static int run_search(long iters, int restarts, unsigned seed) {
    long bound = (long)(n - 1) * (n - 1) / 4;
    srand(seed);
    long best_overall = -1;
    for (int r = 0; r < restarts; r++) {
        for (int i = 0; i < n; i++) pi_[i] = i;
        for (int i = n - 1; i > 0; i--) { int j = rand() % (i + 1); int t = pi_[i]; pi_[i] = pi_[j]; pi_[j] = t; }
        long cur = compute_I();
        double T0 = n / 2.0;
        for (long it = 0; it < iters; it++) {
            double T = T0 * exp(-3.0 * it / iters) + 0.05;
            double u = (double)rand() / RAND_MAX;
            if (u < 0.6) {
                int i = rand() % n, j = rand() % n;
                if (i == j) continue;
                int t = pi_[i]; pi_[i] = pi_[j]; pi_[j] = t;
                long val = compute_I();
                if (val >= cur || (double)rand() / RAND_MAX < exp((val - cur) / T)) cur = val;
                else { t = pi_[i]; pi_[i] = pi_[j]; pi_[j] = t; }
            } else {
                int i = rand() % n, j = rand() % n;
                if (i > j) { int t = i; i = j; j = t; }
                if (i == j) continue;
                for (int l = i, rr = j; l < rr; l++, rr--) { int t = pi_[l]; pi_[l] = pi_[rr]; pi_[rr] = t; }
                long val = compute_I();
                if (val >= cur || (double)rand() / RAND_MAX < exp((val - cur) / T)) cur = val;
                else for (int l = i, rr = j; l < rr; l++, rr--) { int t = pi_[l]; pi_[l] = pi_[rr]; pi_[rr] = t; }
            }
            if (cur > best_overall) best_overall = cur;
        }
    }
    printf("n=%d bound=%ld SA_best=%ld iters=%ld restarts=%d seed=%u %s\n", n, bound, best_overall, iters, restarts,
           seed, best_overall <= bound ? "OK" : "VIOLATION");
    return best_overall <= bound ? 0 : 1;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s exhaustive N | family N [step] | search N iters restarts seed\n", argv[0]);
        return 2;
    }
    n = atoi(argv[2]);
    pi_ = malloc(sizeof(int) * n);
    inv_pi = malloc(sizeof(int) * n);
    a0_vals = malloc(sizeof(int) * n);
    int rc;
    if (strcmp(argv[1], "exhaustive") == 0) {
        rc = run_exhaustive();
    } else if (strcmp(argv[1], "family") == 0) {
        int step = argc >= 4 ? atoi(argv[3]) : 0;
        rc = run_family(step);
    } else if (strcmp(argv[1], "search") == 0) {
        long iters = argc >= 4 ? atol(argv[3]) : 10000;
        int restarts = argc >= 5 ? atoi(argv[4]) : 3;
        unsigned seed = argc >= 6 ? (unsigned)atol(argv[5]) : 1;
        rc = run_search(iters, restarts, seed);
    } else {
        fprintf(stderr, "unknown mode %s\n", argv[1]);
        rc = 2;
    }
    free(pi_); free(inv_pi); free(a0_vals);
    return rc;
}
