/* h13i_value_rotation_only.c — session 9 attempt at H13-I (PLAN §0,
 * docs/notes/h13_line_model.md §6/§7).
 *
 * For every pi of n, restrict the double-cut search to q = 0 (no position
 * cut: read pi in its original index order) and vary only the value shift
 * v0 over all n choices. Report max over all pi of the resulting minimum
 * inversion count (candidate "value_rotation_only" from
 * experiments/h13i_cut_independence.py), against floor((n-1)^2/4).
 *
 * This tests whether the value-rotation freedom alone (without the
 * position-cut freedom) already suffices to meet the H13-I bound. It does
 * not: the excess over the bound is 0 at n=4..6, but 1 at n=7..9 and 3 at
 * n=10 (non-monotone), so this restricted rule is refuted as a route to
 * H13-I. Kept for reproducibility of that negative result.
 *
 * Usage: h13i_value_rotation_only n   (4 <= n <= 11)
 * Version h13i_value_rotation_only-1.0.
 */
#include <stdio.h>
#include <stdlib.h>

static int n;
static int perm[16];
static int used[16];
static int bound_val;
static int max_g;
static int example[16];
static long long count_over_bound;

static int inv_count(const int *seq) {
    int c = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (seq[i] > seq[j]) c++;
    return c;
}

static void process(void) {
    int line[16];
    int best = 1 << 30;
    for (int v0 = 0; v0 < n; v0++) {
        for (int j = 0; j < n; j++) {
            int val = perm[j] - v0;
            if (val < 0) val += n;
            line[j] = val;
        }
        int iv = inv_count(line);
        if (iv < best) best = iv;
    }
    if (best > max_g) {
        max_g = best;
        for (int i = 0; i < n; i++) example[i] = perm[i];
    }
    if (best > bound_val) count_over_bound++;
}

static void rec(int pos) {
    if (pos == n) { process(); return; }
    for (int v = 0; v < n; v++) {
        if (!used[v]) {
            used[v] = 1; perm[pos] = v;
            rec(pos + 1);
            used[v] = 0;
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s n\n", argv[0]); return 2; }
    n = atoi(argv[1]);
    if (n < 2 || n > 11) { fprintf(stderr, "n out of range\n"); return 2; }
    bound_val = ((n - 1) * (n - 1)) / 4;
    max_g = -1;
    count_over_bound = 0;
    rec(0);
    printf("n=%d bound=%d max_over_all_pi=%d excess=%d count_pi_over_bound=%lld example=",
           n, bound_val, max_g, max_g - bound_val, count_over_bound);
    for (int i = 0; i < n; i++) printf("%d ", example[i]);
    printf("\n");
    return 0;
}
