/* h13i_avg_bound.c — exhaustive check of Lemma E (session 9, H13-I attempt):
 *
 *   avg_{a=0}^{n-1} h(a; pi) <= floor((n-1)^2/4)   for every pi in S_n,
 *
 * where, writing pi as a cyclic sequence and beta_a(k) = pi[(k+a+1) mod n]
 * (the "line" obtained by cutting the position circle right after a):
 *
 *   h(a; pi) = min_{b in Z_n} inv( (beta_a(k) - b) mod n )
 *
 * is the minimum inversion count over all VALUE shifts alone (cf.
 * docs/notes/h13_line_model.md, Lemma D below). By definition
 * I(pi) = min_{a,b} inv(...) = min_a h(a; pi), so
 * min_a h(a;pi) <= avg_a h(a;pi); Lemma E, if proved for all n, gives
 * H13-I: I(pi) <= floor((n-1)^2/4).
 *
 * Lemma D (proved, see the note): for fixed a, writing pos = beta_a^{-1},
 *   h(a) = inv(beta_a) - 2 * max_{0<=k<=n} [ P(k) - k(n-1)/2 ],
 *   P(k) = sum_{v=0}^{k-1} pos(v).
 * This is what is computed below (doubled to stay in integers), giving an
 * O(n) per rotation / O(n^2) per permutation algorithm for I(pi) instead of
 * the O(n^4) brute force over all n^2 cuts used by line_profile.c. The two
 * are cross-checked to agree exhaustively for n <= 8 by
 * checks/check_C37_C38.py; this file only checks Lemma E.
 *
 * Usage: h13i_avg_bound n
 *   Exhaustive over all n! permutations (Lehmer order). Reports, in exact
 *   integer arithmetic (everything doubled to avoid the /2 in (n-1)/2):
 *     - the maximum of avg_a h(a) over all pi, compared to floor((n-1)^2/4);
 *     - the number of pi attaining that maximum (expected: exactly n, the
 *       rotations of a single reflection class -- rotating pi permutes the
 *       n values h(a) among themselves, so the sum/avg is rotation-invariant
 *       and reflections come in orbits of size n under rotation of pi itself
 *       composed with the h-in dex rotation).
 *   Also cross-checks Lemma D against an O(n^4) brute force over all (q,c)
 *   cuts (line_profile.c convention) for every permutation, aborting on any
 *   mismatch.
 *
 * Version h13i_avg_bound-1.0.
 */
#include <stdio.h>
#include <stdlib.h>

#define MAXN 11

static int next_perm(int *a, int n) {
    int i = n - 2;
    while (i >= 0 && a[i] >= a[i + 1]) i--;
    if (i < 0) return 0;
    int j = n - 1;
    while (a[j] <= a[i]) j--;
    int t = a[i]; a[i] = a[j]; a[j] = t;
    for (int l = i + 1, r = n - 1; l < r; l++, r--) { t = a[l]; a[l] = a[r]; a[r] = t; }
    return 1;
}

static int inversions(const int *w, int n) {
    int inv = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) if (w[i] > w[j]) inv++;
    return inv;
}

/* Brute force I(pi) over all n^2 cuts, matching line_model.py / line_profile.c. */
static int brute_I(const int *pi, int n) {
    int w[MAXN];
    int best = 1 << 30;
    for (int q = 0; q < n; q++) {
        for (int c = 0; c < n; c++) {
            int shift = q + 1 - c;
            for (int j = 0; j < n; j++) {
                int v = (pi[(q + 1 + j) % n] - shift) % n;
                if (v < 0) v += n;
                w[j] = v;
            }
            int inv = inversions(w, n);
            if (inv < best) best = inv;
        }
    }
    return best;
}

/* 2*h(a) for a single rotation a, and the corresponding beta_a, via Lemma D. */
static long two_h_of_rotation(const int *pi, int n, int a) {
    int beta[MAXN], pos[MAXN];
    for (int k = 0; k < n; k++) beta[k] = pi[(k + a + 1) % n];
    for (int k = 0; k < n; k++) pos[beta[k]] = k;
    int inv = inversions(beta, n);
    long P2 = 0;          /* 2*P(k) */
    long bestD2 = 0;       /* k = 0 term */
    for (int k = 1; k <= n; k++) {
        P2 += 2 * pos[k - 1];
        long D2 = P2 - (long)k * (n - 1);
        if (D2 > bestD2) bestD2 = D2;
    }
    return 2L * inv - 2L * bestD2;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s n\n", argv[0]); return 2; }
    int n = atoi(argv[1]);
    if (n < 2 || n > MAXN) { fprintf(stderr, "n out of range\n"); return 2; }
    int check_lemma_d = (n <= 8); /* brute force O(n^4) per pi -- keep cheap */

    int pi[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;
    long total = 1; for (int k = 2; k <= n; k++) total *= k;

    long bound = (long)(n - 1) * (n - 1) / 4;      /* floor((n-1)^2/4) */
    long best_sum2 = -1;                            /* max over pi of sum_a 2h(a) */
    int best_pi[MAXN];
    long ties = 0;
    long lemma_d_mismatches = 0;
    long checked = 0;

    do {
        if (check_lemma_d) {
            long min_h2 = 1L << 40;
            for (int a = 0; a < n; a++) {
                long h2 = two_h_of_rotation(pi, n, a);
                if (h2 < min_h2) min_h2 = h2;
            }
            int Id = brute_I(pi, n);
            if (min_h2 != 2L * Id) {
                lemma_d_mismatches++;
                fprintf(stderr, "Lemma D mismatch at pi=(");
                for (int i = 0; i < n; i++) fprintf(stderr, "%d%s", pi[i], i + 1 < n ? "," : "");
                fprintf(stderr, "): min_a 2h(a)=%ld  2*brute_I=%d\n", min_h2, 2 * Id);
            }
        }
        long sum2 = 0;
        for (int a = 0; a < n; a++) sum2 += two_h_of_rotation(pi, n, a);
        if (sum2 > best_sum2) {
            best_sum2 = sum2;
            for (int i = 0; i < n; i++) best_pi[i] = pi[i];
            ties = 1;
        } else if (sum2 == best_sum2) {
            ties++;
        }
        checked++;
    } while (next_perm(pi, n));

    double avg_max = (double)best_sum2 / (2.0 * n);
    long rhs2n = 2L * n * bound; /* 2n * bound, compare to best_sum2 */

    printf("{\"n\": %d, \"checked\": %ld, \"lemma_d_checked\": %s, \"lemma_d_mismatches\": %ld,\n",
           n, checked, check_lemma_d ? "true" : "false", lemma_d_mismatches);
    printf(" \"bound_floor_(n-1)^2/4\": %ld, \"max_sum_a_2h(a)\": %ld, \"2n*bound\": %ld,\n",
           bound, best_sum2, rhs2n);
    printf(" \"max_avg_a_h(a)\": %.6f, \"exceeds_bound\": %s, \"num_pi_attaining_max\": %ld,\n",
           avg_max, best_sum2 > rhs2n ? "true" : "false", ties);
    printf(" \"argmax_pi\": [");
    for (int i = 0; i < n; i++) printf("%d%s", best_pi[i], i + 1 < n ? ", " : "");
    printf("]}\n");
    return lemma_d_mismatches ? 1 : 0;
}
