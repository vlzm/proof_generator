/* h13i_balanced_cut-1.0
 *
 * Experiment for H13-I (see docs/notes/h13_line_model.md).
 *
 * Claim (A) "balanced cut exists": for every permutation pi of Z_n there is a
 * double cut (q,s) such that in the resulting line every element takes part in
 * at most floor((n-1)/2) inversions.
 *
 * Claim (B) "counting lemma": every word w of length n in which every element
 * takes part in at most floor((n-1)/2) inversions has
 * inv(w) <= floor((n-1)^2/4).
 *
 * Claim (A') : the inv-minimising cut is itself balanced.
 *
 * Usage: h13i_balanced_cut A <n>   -- claim A and A' for all pi with pi(0)=0
 *        h13i_balanced_cut B <n>   -- claim B for all words of length n
 *
 * Claim A is checked on representatives pi(0)=0 only: the toric class of any
 * permutation contains such a representative (rotate values), and both the
 * claim and inv() are invariant along the class.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int n, thr, bound;

/* per-element inversion counts of the line; returns inv(w), fills r[] */
static int line_stats(const int *w, int *r)
{
    int i, j, inv = 0;
    for (i = 0; i < n; i++) r[i] = 0;
    for (i = 0; i < n; i++)
        for (j = i + 1; j < n; j++)
            if (w[i] > w[j]) { r[i]++; r[j]++; inv++; }
    return inv;
}

static long long count_perm;
static long long count_bad_A, count_bad_Aprime, count_bad_B;
static int worst_maxinv_A;      /* max over pi of min over cuts of max_v inv_v */
static int worst_inv_B;         /* max inv over balanced words */
static int example_shown;

static void check_A(const int *pi)
{
    int q, s, i, w[16], r[16];
    int best_max = 1 << 30, best_inv = 1 << 30, best_inv_max = 0;
    count_perm++;
    for (q = 0; q < n; q++)
        for (s = 0; s < n; s++) {
            int inv, mx = 0;
            for (i = 0; i < n; i++) {
                int v = pi[(q + i) % n] - s;
                w[i] = v < 0 ? v + n : v;
            }
            inv = line_stats(w, r);
            for (i = 0; i < n; i++) if (r[i] > mx) mx = r[i];
            if (mx < best_max) best_max = mx;
            if (inv < best_inv) { best_inv = inv; best_inv_max = mx; }
            else if (inv == best_inv && mx < best_inv_max) best_inv_max = mx;
        }
    if (best_max > thr) {
        count_bad_A++;
        if (example_shown < 5) {
            printf("  A FAILS: pi =");
            for (i = 0; i < n; i++) printf(" %d", pi[i]);
            printf("  min_cut max_v inv_v = %d > %d\n", best_max, thr);
            example_shown++;
        }
    }
    if (best_max > worst_maxinv_A) worst_maxinv_A = best_max;
    if (best_inv_max > thr) {
        count_bad_Aprime++;
        if (example_shown < 10) {
            printf("  A' FAILS: pi =");
            for (i = 0; i < n; i++) printf(" %d", pi[i]);
            printf("  inv-min cut has max_v inv_v = %d > %d\n", best_inv_max, thr);
            example_shown++;
        }
    }
}

static void check_B(const int *w)
{
    int r[16], i, mx = 0, inv;
    count_perm++;
    inv = line_stats(w, r);
    for (i = 0; i < n; i++) if (r[i] > mx) mx = r[i];
    if (mx > thr) return;
    if (inv > worst_inv_B) worst_inv_B = inv;
    if (inv > bound) {
        count_bad_B++;
        if (example_shown < 5) {
            printf("  B FAILS: w =");
            for (i = 0; i < n; i++) printf(" %d", w[i]);
            printf("  inv = %d > %d\n", inv, bound);
            example_shown++;
        }
    }
}

/* ---- mode S: max of S over the n^2 cuts vs max over all 2^(n-1) switchings */

static long long gap_hist[64];
static int gap_max;
static long long count_S_bad;   /* max over cuts < floor(n/2) : would refute H13-I */

static void check_S(const int *pi)
{
    int chi[16][16];            /* signing at the cut (0,0) */
    int q, s, i, j, u, v, w[16];
    int best_cut = -1000, best_all = -1000, N = n * (n - 1) / 2;
    count_perm++;
    /* base signing: element index i -> (line position, line value) at cut (0,0) */
    {
        int posn[16], val[16];
        for (i = 0; i < n; i++) { posn[i] = i; val[i] = pi[i]; }
        for (u = 0; u < n; u++)
            for (v = 0; v < n; v++)
                if (u != v)
                    chi[u][v] = ((posn[u] - posn[v]) * (val[u] - val[v]) > 0) ? 1 : -1;
    }
    for (q = 0; q < n; q++)
        for (s = 0; s < n; s++) {
            int inv, S, r[16];
            for (i = 0; i < n; i++) {
                int x = pi[(q + i) % n] - s;
                w[i] = x < 0 ? x + n : x;
            }
            inv = line_stats(w, r);
            S = N - 2 * inv;
            if (S > best_cut) best_cut = S;
        }
    for (j = 0; j < (1 << (n - 1)); j++) {
        int S = 0;
        for (u = 0; u < n; u++)
            for (v = u + 1; v < n; v++) {
                int mu = (u < n - 1) ? (j >> u) & 1 : 0;
                int mv = (v < n - 1) ? (j >> v) & 1 : 0;
                S += (mu != mv) ? -chi[u][v] : chi[u][v];
            }
        if (S > best_all) best_all = S;
    }
    if (best_cut < n / 2) {
        count_S_bad++;
        if (example_shown < 5) {
            printf("  H13-I FAILS: pi =");
            for (i = 0; i < n; i++) printf(" %d", pi[i]);
            printf("  max_cut S = %d < %d\n", best_cut, n / 2);
            example_shown++;
        }
    }
    {
        int g = best_all - best_cut;
        if (g > gap_max) gap_max = g;
        if (g < 64) gap_hist[g]++;
        if (g > 0 && example_shown < 8) {
            printf("  gap %d: pi =", g);
            for (i = 0; i < n; i++) printf(" %d", pi[i]);
            printf("  max_cut S = %d, max_switch S = %d\n", best_cut, best_all);
            example_shown++;
        }
    }
}

/* enumerate permutations of 0..n-1 (optionally with a[0]=0 fixed) */
static void rec(int *a, int *used, int k, int fix0, void (*cb)(const int *))
{
    int v;
    if (k == n) { cb(a); return; }
    for (v = 0; v < n; v++) {
        if (used[v]) continue;
        if (k == 0 && fix0 && v != 0) continue;
        used[v] = 1; a[k] = v;
        rec(a, used, k + 1, fix0, cb);
        used[v] = 0;
    }
}

int main(int argc, char **argv)
{
    int a[16], used[16];
    char mode;
    if (argc < 3) { fprintf(stderr, "usage: %s A|B n\n", argv[0]); return 2; }
    mode = argv[1][0];
    n = atoi(argv[2]);
    if (n < 2 || n > 13) { fprintf(stderr, "n out of range\n"); return 2; }
    thr = (n - 1) / 2;
    bound = (n - 1) * (n - 1) / 4;
    memset(used, 0, sizeof used);
    if (mode == 'A') {
        rec(a, used, 0, 1, check_A);
        printf("A n=%d thr=%d perms(pi0=0)=%lld bad_A=%lld bad_Aprime=%lld "
               "max_pi min_cut max_v inv_v=%d\n",
               n, thr, count_perm, count_bad_A, count_bad_Aprime, worst_maxinv_A);
    } else if (mode == 'S') {
        int g;
        rec(a, used, 0, 1, check_S);
        printf("S n=%d perms(pi0=0)=%lld H13I_violations=%lld max_gap=%d gaps:",
               n, count_perm, count_S_bad, gap_max);
        for (g = 0; g <= gap_max && g < 64; g++)
            printf(" %d:%lld", g, gap_hist[g]);
        printf("\n");
    } else {
        rec(a, used, 0, 0, check_B);
        printf("B n=%d thr=%d bound=%d words=%lld bad_B=%lld max_inv_balanced=%d\n",
               n, thr, bound, count_perm, count_bad_B, worst_inv_B);
    }
    return 0;
}
