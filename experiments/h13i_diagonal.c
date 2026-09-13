/* h13i_diagonal.c -- exhaustive tests of three candidate reductions of H13-I
 * (I(pi) <= floor((n-1)^2/4) for every pi, via the double cut (q, c)) to a
 * single free parameter, plus a local-move (transposition) test.
 *
 * For every permutation pi of 0..n-1 (Lehmer order), with M = floor((n-1)^2/4):
 *
 *  (1) I_shift0(pi)  = min over q of inv(w) on the diagonal s = 0 (c = q+1),
 *      i.e. values unchanged, only the starting point of the linear read
 *      varies (Lemma D of docs/notes/h13_line_model.md #7).
 *  (2) I_fixedc(pi)  = min over c of inv(w) with q fixed at n-1 (natural
 *      position order, only the value labelling is rotated) -- the dual of
 *      (1) under pi -> pi^{-1}.
 *  (3) minqA(pi)     = min over q of the exact average over c of inv(w)
 *      (n * minqA is an integer, reported as such): does choosing q to
 *      minimise the c-average already reach the bound?
 *  (4) I_full(pi)    = min over all n^2 cuts (q, c) -- ground truth, used to
 *      flag permutations with I_full < M that admit NO transposition (swap
 *      of two values pi(i), pi(j)) strictly increasing I_full ("stuck" local
 *      maxima under the transposition move; only computed for n <= LOCAL_NMAX
 *      because of the O(n^2) transpositions x O(n^2) cuts x O(n^2) inversions
 *      cost per permutation).
 *
 * Usage: h13i_diagonal n [--local]   (--local enables test (4); n <= LOCAL_NMAX)
 * Output: one JSON object on stdout. Version h13i_diagonal-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 12
#define LOCAL_NMAX 7

static int inversions(const int *w, int n) {
    int inv = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (w[i] > w[j]) inv++;
    return inv;
}

static int line_inv(const int *pi, int n, int q, int c) {
    int w[MAXN];
    int shift = q + 1 - c;
    for (int j = 0; j < n; j++) {
        int v = (pi[(q + 1 + j) % n] - shift) % n;
        if (v < 0) v += n;
        w[j] = v;
    }
    return inversions(w, n);
}

static int I_full(const int *pi, int n) {
    int best = 999999;
    for (int q = 0; q < n; q++)
        for (int c = 0; c < n; c++) {
            int iv = line_inv(pi, n, q, c);
            if (iv < best) best = iv;
        }
    return best;
}

static int I_shift0(const int *pi, int n) {
    int best = 999999;
    for (int q = 0; q < n; q++) {
        int c = (q + 1) % n; /* shift = q + 1 - c = 0 */
        int iv = line_inv(pi, n, q, c);
        if (iv < best) best = iv;
    }
    return best;
}

static int I_fixedc(const int *pi, int n) {
    int best = 999999;
    int q = n - 1;
    for (int c = 0; c < n; c++) {
        int iv = line_inv(pi, n, q, c);
        if (iv < best) best = iv;
    }
    return best;
}

/* n * (min over q of average over c of inv(w)) -- exact integer */
static long minqA_times_n(const int *pi, int n) {
    long best = -1;
    for (int q = 0; q < n; q++) {
        long sum = 0;
        for (int c = 0; c < n; c++) sum += line_inv(pi, n, q, c);
        if (best < 0 || sum < best) best = sum;
    }
    return best;
}

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

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s n [--local]\n", argv[0]); return 2; }
    int n = atoi(argv[1]);
    int do_local = (argc >= 3 && strcmp(argv[2], "--local") == 0);
    if (n < 2 || n > MAXN) return 2;
    if (do_local && n > LOCAL_NMAX) { fprintf(stderr, "--local only for n <= %d\n", LOCAL_NMAX); return 2; }

    long M = (long)(n - 1) * (n - 1) / 4;

    int pi[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;

    long total = 1;
    for (int k = 2; k <= n; k++) total *= k;

    int max_shift0 = -1, max_fixedc = -1;
    long argmax_shift0 = -1, argmax_fixedc = -1;
    long max_minqA_num = -1, argmax_minqA = -1; /* numerator (times n) */
    long count_shift0_exceeds = 0, count_fixedc_exceeds = 0, count_minqA_exceeds = 0;
    long count_stuck = 0; /* I_full < M with no improving transposition */
    int stuck_examples[5][MAXN]; int stuck_example_I[5]; int n_stuck_examples = 0;

    long rank = 0;
    do {
        int s0 = I_shift0(pi, n);
        int fc = I_fixedc(pi, n);
        long mA = minqA_times_n(pi, n);
        if (s0 > max_shift0) { max_shift0 = s0; argmax_shift0 = rank; }
        if (fc > max_fixedc) { max_fixedc = fc; argmax_fixedc = rank; }
        if (mA > max_minqA_num) { max_minqA_num = mA; argmax_minqA = rank; }
        if (s0 > M) count_shift0_exceeds++;
        if (fc > M) count_fixedc_exceeds++;
        if (mA > M * (long)n) count_minqA_exceeds++;

        if (do_local) {
            int ifull = I_full(pi, n);
            if (ifull < M) {
                int improved = 0;
                int p[MAXN];
                memcpy(p, pi, n * sizeof(int));
                for (int i = 0; i < n && !improved; i++)
                    for (int j = i + 1; j < n && !improved; j++) {
                        int t = p[i]; p[i] = p[j]; p[j] = t;
                        if (I_full(p, n) > ifull) improved = 1;
                        t = p[i]; p[i] = p[j]; p[j] = t;
                    }
                if (!improved) {
                    count_stuck++;
                    if (n_stuck_examples < 5) {
                        memcpy(stuck_examples[n_stuck_examples], pi, n * sizeof(int));
                        stuck_example_I[n_stuck_examples] = ifull;
                        n_stuck_examples++;
                    }
                }
            }
        }
        rank++;
    } while (next_perm(pi, n));

    printf("{\"n\": %d, \"count\": %ld, \"M\": %ld,\n", n, total, M);
    printf(" \"shift0\": {\"max\": %d, \"argmax_rank\": %ld, \"count_exceeds_M\": %ld},\n",
           max_shift0, argmax_shift0, count_shift0_exceeds);
    printf(" \"fixedc\": {\"max\": %d, \"argmax_rank\": %ld, \"count_exceeds_M\": %ld},\n",
           max_fixedc, argmax_fixedc, count_fixedc_exceeds);
    printf(" \"minqA\": {\"max_numerator_times_n\": %ld, \"max_value\": %.6f, \"argmax_rank\": %ld, \"count_exceeds_M\": %ld}",
           max_minqA_num, (double)max_minqA_num / n, argmax_minqA, count_minqA_exceeds);
    if (do_local) {
        printf(",\n \"local_transposition\": {\"count_stuck_below_M\": %ld, \"examples\": [", count_stuck);
        for (int e = 0; e < n_stuck_examples; e++) {
            printf("%s{\"pi\": [", e ? ", " : "");
            for (int i = 0; i < n; i++) printf("%s%d", i ? "," : "", stuck_examples[e][i]);
            printf("], \"I\": %d}", stuck_example_I[e]);
        }
        printf("]}");
    }
    printf("\n}\n");
    return 0;
}
