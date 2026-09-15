/* h13i_negative_results.c — session 9: two more rejected approaches to H13-I
 * (docs/notes/h13_line_model.md #6: "for every pi there is a double cut
 * (position rotation q, value rotation c) with inv(w) <= floor((n-1)^2/4)").
 *
 * I(pi) = min over q, c in Z_n of inv of the relabelled line (see
 * experiments/line_model.py, experiments/line_profile.c). This program checks,
 * exhaustively over all permutations of n, two weaker claims that would have
 * made H13-I easy if true:
 *
 * (A) "Domain rotation alone suffices": fix q = n-1 (read pi in its original
 *     order) and minimize only over the value rotation c.  Tests whether
 *     max_pi [ (min_c inv) - floor((n-1)^2/4) ] <= 0.
 *
 * (B) "Pair-removal induction": for every pi of size n there is a pair of
 *     positions {p, r} such that contracting them (delete positions p, r and
 *     values pi(p), pi(r); relabel remaining positions and values 0..n-3
 *     preserving order) gives a size-(n-2) permutation pi' with
 *       I(pi) - I(pi') <= floor((n-1)^2/4) - floor((n-3)^2/4).
 *     This is the natural 2-element generalisation of the single-element
 *     induction already rejected in h13_line_model.md #1.7.  Tests whether
 *     max_pi [ min_{p,r} (I(pi) - I(pi')) - budget(n) ] <= 0, and how the gap
 *     behaves as n grows (single-element induction failed at one n = 8
 *     example without checking growth; this program also reports the worst
 *     gap so growth can be read off across n).
 *
 * Both (A) and (B) are found FALSE already at small n (n = 6 and n = 8
 * respectively), so neither collapses the n^2 search.  See also
 * experiments/h13i_affine_growth.py for the same test (B) on the affine
 * family pi(i) = a*i + b mod n at larger n, where the gap grows with n
 * rather than being a small-n artefact.
 *
 * Usage: gcc -O2 -o h13i_negative_results experiments/h13i_negative_results.c
 *        ./h13i_negative_results nmax
 * Exhaustive over all n! permutations for 4 <= n <= nmax (Lehmer/lexicographic
 * order).  No dependency on oracle/moves.py: I(pi) is a purely combinatorial
 * quantity (already defined and used without the oracle in line_profile.c).
 * Version h13i_negative_results-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static long long fact(int k) { long long r = 1; for (int i = 2; i <= k; i++) r *= i; return r; }

static long long rank_perm(const int *p, int n) {
    long long rank = 0;
    int used[16] = {0};
    for (int i = 0; i < n; i++) {
        int cnt = 0;
        for (int v = 0; v < p[i]; v++) if (!used[v]) cnt++;
        rank += (long long)cnt * fact(n - 1 - i);
        used[p[i]] = 1;
    }
    return rank;
}

static void unrank_perm(long long rank, int n, int *out) {
    int avail[16];
    for (int i = 0; i < n; i++) avail[i] = i;
    int navail = n;
    for (int i = 0; i < n; i++) {
        long long f = fact(n - 1 - i);
        int idx = (int)(rank / f);
        rank %= f;
        out[i] = avail[idx];
        for (int j = idx; j < navail - 1; j++) avail[j] = avail[j + 1];
        navail--;
    }
}

static int inv_count(const int *a, int n) {
    int c = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (a[i] > a[j]) c++;
    return c;
}

/* I(pi): min over q (domain rotation) and c (value rotation) of inv of the
 * relabelled line; equivalent to min over independent domain shift s and
 * value shift t of inv({(pi[(s+i)%n]-t) mod n}_i) -- see line_model.py. */
static int I_of(const int *pi, int n) {
    int best = 1 << 30;
    int rot[16], w[16];
    for (int s = 0; s < n; s++) {
        for (int i = 0; i < n; i++) rot[i] = pi[(s + i) % n];
        for (int t = 0; t < n; t++) {
            for (int i = 0; i < n; i++) w[i] = ((rot[i] - t) % n + n) % n;
            int iv = inv_count(w, n);
            if (iv < best) best = iv;
        }
    }
    return best;
}

/* (A): min over c only, s fixed to 0 (read pi in its given order). */
static int I_fixed_domain(const int *pi, int n) {
    int best = 1 << 30;
    int w[16];
    for (int t = 0; t < n; t++) {
        for (int i = 0; i < n; i++) w[i] = ((pi[i] - t) % n + n) % n;
        int iv = inv_count(w, n);
        if (iv < best) best = iv;
    }
    return best;
}

static void contract(const int *pi, int n, int p, int r, int *out) {
    int vals[16], nv = 0;
    for (int i = 0; i < n; i++) if (i != p && i != r) vals[nv++] = pi[i];
    int sorted[16];
    memcpy(sorted, vals, nv * sizeof(int));
    for (int i = 1; i < nv; i++) {
        int key = sorted[i], j = i - 1;
        while (j >= 0 && sorted[j] > key) { sorted[j + 1] = sorted[j]; j--; }
        sorted[j + 1] = key;
    }
    for (int i = 0; i < nv; i++) {
        int rk = 0;
        for (int k = 0; k < nv; k++) if (sorted[k] < vals[i]) rk++;
        out[i] = rk;
    }
}

int main(int argc, char **argv) {
    int nmax = argc > 1 ? atoi(argv[1]) : 9;
    printf("{\"version\": \"h13i_negative_results-1.0\", \"rows\": [\n");
    int first_row = 1;
    for (int n = 4; n <= nmax; n++) {
        long long total = fact(n);
        int bound_n = ((n - 1) * (n - 1)) / 4;

        /* precompute I(pi') table for size m = n - 2, if n >= 6 */
        int m = n - 2;
        unsigned char *Itab = NULL;
        int bound_m = 0, budget = 0;
        if (m >= 4) {
            long long totalm = fact(m);
            Itab = malloc(totalm);
            int perm[16];
            for (long long r = 0; r < totalm; r++) {
                unrank_perm(r, m, perm);
                Itab[r] = (unsigned char)I_of(perm, m);
            }
            bound_m = ((m - 1) * (m - 1)) / 4;
            budget = bound_n - bound_m;
        }

        int maxA_gap = -1; long long argA = -1;
        long long failA = 0;
        int maxB_needed = -1; long long argB = -1;
        long long failB = 0;
        int perm[16], reduced[16];

        for (long long r = 0; r < total; r++) {
            unrank_perm(r, n, perm);

            int Ifixed = I_fixed_domain(perm, n);
            int gapA = Ifixed - bound_n;
            if (gapA > maxA_gap) { maxA_gap = gapA; argA = r; }
            if (gapA > 0) failA++;

            if (Itab) {
                int Ipi = I_of(perm, n);
                int best_reduced = 1 << 30;
                for (int p = 0; p < n; p++)
                    for (int rr = p + 1; rr < n; rr++) {
                        contract(perm, n, p, rr, reduced);
                        long long rk = rank_perm(reduced, m);
                        if (Itab[rk] < best_reduced) best_reduced = Itab[rk];
                    }
                int needed = Ipi - best_reduced;
                if (needed > maxB_needed) { maxB_needed = needed; argB = r; }
                if (needed > budget) failB++;
            }
        }

        int argA_perm[16], argB_perm[16];
        unrank_perm(argA, n, argA_perm);
        if (!first_row) printf(",\n");
        first_row = 0;
        printf("  {\"n\": %d, \"bound_n\": %d, \"checked\": %lld,\n", n, bound_n, total);
        printf("   \"testA_fixed_domain_max_gap\": %d, \"testA_fails\": %lld, \"testA_argmax\": [",
               maxA_gap, failA);
        for (int i = 0; i < n; i++) printf("%d%s", argA_perm[i], i + 1 < n ? ", " : "");
        printf("],\n");
        if (Itab) {
            unrank_perm(argB, n, argB_perm);
            printf("   \"testB_pair_removal_budget\": %d, \"testB_max_needed\": %d, \"testB_fails\": %lld, \"testB_argmax\": [",
                   budget, maxB_needed, failB);
            for (int i = 0; i < n; i++) printf("%d%s", argB_perm[i], i + 1 < n ? ", " : "");
            printf("]\n");
            free(Itab);
        } else {
            printf("   \"testB_pair_removal_budget\": null, \"testB_max_needed\": null, \"testB_fails\": null, \"testB_argmax\": null\n");
        }
        printf("  }");
        fflush(stdout);
    }
    printf("\n]}\n");
    return 0;
}
