/* toric_local_conditions.c — how deep must the local optimality conditions of a
 * double cut be before they force inv(w) <= floor((n-1)^2/4)?
 *
 * For a line w (a permutation of Z_n) the shift of the cut by (k, l) changes the
 * inversion count by (docs/proofs/C37_toric_inversions.md, Cor. 1.7)
 *
 *   Delta(k, l) = A(k) + B(l) + 2 k l - 2 n K(k, l),
 *   A(k) = sum_{j<k} (n - 1 - 2 w_j),      B(l) = sum_{y<l} (n - 1 - 2 t_y),
 *   K(k, l) = #{j < k : w_j < l},          t_y = position of the value y,
 *
 * so "the cut is optimal in its toric class" means Delta(k, l) >= 0 for all
 * 0 <= k, l <= n.  The program enumerates all lines that satisfy the subfamily
 * of these conditions with min(k, l) <= M ("depth M"; M = 0 is Lemma C of
 * docs/notes/h13_line_model.md) and reports the maximum of inv(w) over them.
 * Depth 1 suffices for 5 <= n <= 9 and fails for n = 10; depth 2 suffices for
 * 5 <= n <= 10.  If a fixed depth sufficed for all n, H13-I would follow from
 * a finite inequality; the n = 10 counterexample shows depth 1 is not enough.
 *
 * Usage: gcc -O2 -o toric_local_conditions experiments/toric_local_conditions.c
 *        ./toric_local_conditions n M
 * Version toric_local-1.0.
 */
#include <stdio.h>
#include <stdlib.h>

static int n, M;

static int next_perm(int *a, int m) {
    int i = m - 2; while (i >= 0 && a[i] >= a[i + 1]) i--;
    if (i < 0) return 0;
    int j = m - 1; while (a[j] <= a[i]) j--;
    int t = a[i]; a[i] = a[j]; a[j] = t;
    for (int l = i + 1, r = m - 1; l < r; l++, r--) { t = a[l]; a[l] = a[r]; a[r] = t; }
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s n M\n", argv[0]); return 2; }
    n = atoi(argv[1]); M = atoi(argv[2]);
    if (n < 3 || n > 12) return 2;
    int w[16], t[16], A[17], B[17], bw[16];
    for (int i = 0; i < n; i++) w[i] = i;
    int best = -1; long cnt_ok = 0, total = 0;
    do {
        total++;
        for (int j = 0; j < n; j++) t[w[j]] = j;
        A[0] = 0; for (int k = 0; k < n; k++) A[k + 1] = A[k] + n - 1 - 2 * w[k];
        B[0] = 0; for (int l = 0; l < n; l++) B[l + 1] = B[l] + n - 1 - 2 * t[l];
        int ok = 1;
        for (int k = 0; k <= n && ok; k++)
            for (int l = 0; l <= n; l++) {
                if ((k < l ? k : l) > M) continue;
                int K = 0; for (int j = 0; j < k; j++) if (w[j] < l) K++;
                if (A[k] + B[l] + 2 * k * l - 2 * n * K < 0) { ok = 0; break; }
            }
        if (!ok) continue;
        cnt_ok++;
        int inv = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) if (w[i] > w[j]) inv++;
        if (inv > best) { best = inv; for (int i = 0; i < n; i++) bw[i] = w[i]; }
    } while (next_perm(w, n));
    printf("{\"mode\": \"local\", \"n\": %d, \"depth\": %d, \"target\": %d, \"lines\": %ld, "
           "\"lines_satisfying\": %ld, \"max_inv\": %d, \"sufficient\": %s, \"argmax\": [",
           n, M, ((n - 1) * (n - 1)) / 4, total, cnt_ok, best,
           best <= ((n - 1) * (n - 1)) / 4 ? "true" : "false");
    for (int i = 0; i < n; i++) printf("%s%d", i ? ", " : "", bw[i]);
    printf("]}\n");
    return 0;
}
