/* h13i_axis_scan.c — H13-I (PLAN §8): does every pi admit a double cut
 * (edge {q, q+1} on positions, shift c on values) with at most
 * floor((n-1)^2/4) inversions of the relabelled line?  See
 * docs/notes/h13_line_model.md §6 for the statement and rejected approaches.
 *
 * This tool checks three things exhaustively over all pi in S_n:
 *
 *   mode=full   I(pi) = min over ALL n^2 cuts (q, c) of inv(line).
 *               Extends the C33/C35 exhaustive range (previously 4<=n<=10,
 *               experiments/line_profile.c) using the symmetry pi(0) = 0:
 *               I(pi) is invariant under replacing pi by a value-shift of
 *               pi (shifting is exactly part of the c-search), so it
 *               suffices to enumerate the (n-1)! permutations with
 *               pi(0) = 0 and still take max over all of them.
 *
 *   mode=axis   Only the value-shift c is free; the position cut is fixed
 *               at a = 0 (pi read in its given order).  Tests whether
 *               value-rotation alone (no position-rotation) already
 *               reaches floor(n^2/4) non-inversions, i.e. whether the
 *               position-cut freedom in H13-I is dispensable.
 *
 *   mode=diag   Only the n "diagonal" cuts with a = c are tried (position
 *               cut and value shift forced equal), instead of the full n^2.
 *
 * mode=full reports max I over the pi(0)=0 representatives against
 * floor((n-1)^2/4); mode=axis/diag report min-over-pi of (max non-inversions
 * within the restricted cut family) against floor(n^2/4) (a lower target
 * means the restricted family is insufficient) plus the worst pi found.
 *
 * Usage: h13i_axis_scan n mode   (mode in {full, axis, diag})
 * Version h13i_axis_scan-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 16

static int n;
static int perm[MAXN], used[MAXN];
static long long checked;

static inline int inv_count(const int *seq) {
    int c = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (seq[i] > seq[j]) c++;
    return c;
}

/* full double-cut minimum inversions over all n^2 cuts */
static int min_inv_full(const int *pi) {
    int best = 1 << 30;
    int w[MAXN], wp[MAXN];
    for (int a = 0; a < n; a++) {
        for (int k = 0; k < n; k++) w[k] = pi[(a + k) % n];
        for (int b = 0; b < n; b++) {
            for (int k = 0; k < n; k++) wp[k] = (w[k] - b + 100 * n) % n;
            int iv = inv_count(wp);
            if (iv < best) best = iv;
            if (best == 0) return 0;
        }
    }
    return best;
}

/* value-shift only (position order fixed as given): max non-inversions over
 * the n shifts b */
static int max_noninv_axis(const int *pi) {
    int total = n * (n - 1) / 2;
    int best = -1;
    int wp[MAXN];
    for (int b = 0; b < n; b++) {
        for (int k = 0; k < n; k++) wp[k] = (pi[k] - b + 100 * n) % n;
        int noninv = total - inv_count(wp);
        if (noninv > best) best = noninv;
    }
    return best;
}

/* diagonal family a = c: max non-inversions over the n matched cuts */
static int max_noninv_diag(const int *pi) {
    int total = n * (n - 1) / 2;
    int best = -1;
    int w[MAXN], wp[MAXN];
    for (int a = 0; a < n; a++) {
        for (int k = 0; k < n; k++) w[k] = pi[(a + k) % n];
        for (int k = 0; k < n; k++) wp[k] = (w[k] - a + 100 * n) % n;
        int noninv = total - inv_count(wp);
        if (noninv > best) best = noninv;
    }
    return best;
}

/* ---- mode=full: max over pi(0)=0 representatives ---- */
static int full_best;
static long long full_achievers;

static void recurse_full(int depth) {
    if (depth == n) {
        int I = min_inv_full(perm);
        checked++;
        if (I > full_best) { full_best = I; full_achievers = 1; }
        else if (I == full_best) full_achievers++;
        return;
    }
    for (int v = 0; v < n; v++) {
        if (!used[v]) {
            used[v] = 1;
            perm[depth] = v;
            recurse_full(depth + 1);
            used[v] = 0;
        }
    }
}

/* ---- mode=axis / mode=diag: min over ALL pi (no symmetry reduction,
 * since axis/diag are not invariant under an arbitrary value-shift of pi
 * the way the full search is) ---- */
static int restricted_worst;
static long long restricted_worst_count;
static const char *restricted_mode;

static void recurse_restricted(int depth) {
    if (depth == n) {
        int m = (strcmp(restricted_mode, "axis") == 0) ? max_noninv_axis(perm)
                                                          : max_noninv_diag(perm);
        checked++;
        if (m < restricted_worst) { restricted_worst = m; restricted_worst_count = 1; }
        else if (m == restricted_worst) restricted_worst_count++;
        return;
    }
    for (int v = 0; v < n; v++) {
        if (!used[v]) {
            used[v] = 1;
            perm[depth] = v;
            recurse_restricted(depth + 1);
            used[v] = 0;
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s n {full|axis|diag}\n", argv[0]);
        return 1;
    }
    n = atoi(argv[1]);
    if (n < 1 || n > MAXN) { fprintf(stderr, "n out of range\n"); return 1; }
    const char *mode = argv[2];
    memset(used, 0, sizeof(used));
    checked = 0;

    if (strcmp(mode, "full") == 0) {
        full_best = -1;
        full_achievers = 0;
        used[0] = 1;
        perm[0] = 0;
        recurse_full(1);
        int target = (n - 1) * (n - 1) / 4;
        printf("n=%d mode=full checked=%lld max_I=%d target_floor((n-1)^2/4)=%d "
               "achievers=%lld %s\n",
               n, checked, full_best, target, full_achievers,
               full_best <= target ? "OK" : "FAIL");
    } else if (strcmp(mode, "axis") == 0 || strcmp(mode, "diag") == 0) {
        restricted_mode = mode;
        restricted_worst = 1 << 30;
        restricted_worst_count = 0;
        recurse_restricted(0);
        int target = n * n / 4;
        printf("n=%d mode=%s checked=%lld min_over_pi_max_over_cuts_noninv=%d "
               "target_floor(n^2/4)=%d gap=%d worst_count=%lld %s\n",
               n, mode, checked, restricted_worst, target,
               target - restricted_worst, restricted_worst_count,
               restricted_worst >= target ? "OK" : "FAIL");
    } else {
        fprintf(stderr, "unknown mode %s\n", mode);
        return 1;
    }
    return 0;
}
