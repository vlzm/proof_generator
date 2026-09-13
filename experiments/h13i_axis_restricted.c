/* h13i_axis_restricted.c — exhaustive check for H13-I (session 9):
 * for every permutation pi of {0,...,n-1}, is there a double cut (q, c)
 * (position-edge rotation a and value shift b, see experiments/line_model.py
 * and docs/notes/h13_line_model.md §0) with inv(w) <= floor((n-1)^2/4)?
 * Reformulation used here (docs/notes/h13i_attempt.md §1): the double cut
 * is exactly the group action of Z_n x Z_n on permutations of Z_n,
 * pi_{a,b}(j) = (pi((j+a) mod n) - b) mod n, and I(pi) = min_{a,b} inv(pi_{a,b}).
 *
 * This program checks, exhaustively over all pi (all n! of them, generated
 * in lexicographic order), three quantities against bound = floor((n-1)^2/4):
 *   - full search over all (a,b) in Z_n x Z_n (sanity re-check of C33/C35's
 *     already-certified claim I(pi) <= bound; must give 0 failures);
 *   - "b-only": a fixed at 0 (no position rotation), vary only the value
 *     shift b;
 *   - "a-only": b fixed at 0 (no value shift), vary only the position
 *     rotation a (a plain cyclic rotation of the sequence).
 * b-only and a-only are the two natural single-axis restrictions of the
 * double cut; this checks whether either alone would already suffice
 * (it would make H13-I trivial). Independent of the LRX oracle: this is a
 * pure statement about inversions of integer sequences.
 *
 * Usage: h13i_axis_restricted n
 * Version h13i_axis_restricted-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int n;
static int perm[16], used[16];
static int bound_lo;
static long total = 0, b_only_fail = 0, a_only_fail = 0, full_fail = 0;
static int rotv[16], shiftv[16];

static int inversions(const int *w, int m) {
    int inv = 0;
    for (int i = 0; i < m; i++)
        for (int j = i + 1; j < m; j++)
            if (w[i] > w[j]) inv++;
    return inv;
}

static void process(void) {
    total++;

    int best_b = n * n;
    for (int b = 0; b < n; b++) {
        for (int j = 0; j < n; j++) shiftv[j] = ((perm[j] - b) % n + n) % n;
        int iv = inversions(shiftv, n);
        if (iv < best_b) best_b = iv;
    }
    if (best_b > bound_lo) b_only_fail++;

    int best_a = n * n;
    for (int a = 0; a < n; a++) {
        for (int j = 0; j < n; j++) rotv[j] = perm[(j + a) % n];
        int iv = inversions(rotv, n);
        if (iv < best_a) best_a = iv;
    }
    if (best_a > bound_lo) a_only_fail++;

    int best_full = n * n;
    for (int a = 0; a < n; a++) {
        for (int j = 0; j < n; j++) rotv[j] = perm[(j + a) % n];
        for (int b = 0; b < n; b++) {
            for (int j = 0; j < n; j++) shiftv[j] = ((rotv[j] - b) % n + n) % n;
            int iv = inversions(shiftv, n);
            if (iv < best_full) best_full = iv;
        }
    }
    if (best_full > bound_lo) full_fail++;
}

static void gen(int pos) {
    if (pos == n) { process(); return; }
    for (int v = 0; v < n; v++) {
        if (!used[v]) {
            used[v] = 1; perm[pos] = v;
            gen(pos + 1);
            used[v] = 0;
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s n\n", argv[0]); return 2; }
    n = atoi(argv[1]);
    if (n < 2 || n > 11) { fprintf(stderr, "n out of range\n"); return 2; }
    bound_lo = ((n - 1) * (n - 1)) / 4;
    memset(used, 0, sizeof(used));
    gen(0);
    printf("{\"n\": %d, \"bound\": %d, \"total\": %ld, "
           "\"b_only_fail\": %ld, \"a_only_fail\": %ld, \"full_fail\": %ld}\n",
           n, bound_lo, total, b_only_fail, a_only_fail, full_fail);
    return 0;
}
