/* h3_reduction.c -- test the H3 induction hypothesis (PLAN.md H3 / S4.3.3):
 * is there, for every pi in S_n, a way to remove one element (closing the
 * gap in the array and relabelling the remaining n-1 values down) so that
 * the resulting pi' in S_{n-1} satisfies d_n(pi) <= d_{n-1}(pi') + (n-1)?
 * That exact recurrence, with base case d(id_4)=... D_4 = B_4 = 6, would
 * telescope to D_n <= B_n with zero slack (B_n - B_{n-1} = n-1 always).
 *
 * "Remove element v" (v = 0..n-1, n choices) means: delete the array entry
 * equal to v, keep the relative order of the rest, then relabel values
 * greater than v down by one so the result is a genuine permutation of
 * 0..n-2. This is the most natural embedding of S_{n-1} as a "one element
 * fixed" subfamily of S_n; PLAN.md's own warning ("a fixed element cannot
 * be silently dropped -- L and R rotate the whole array") suggests the
 * *position* where the removed array closes up matters, so besides the
 * "natural" head-relative reading order after deletion, this program also
 * tries EVERY rotation of the resulting (n-1)-array (an L/R walk in the
 * (n-1)-model) and reports the minimum distance reachable that way -- this
 * is a strict upper bound on what any single-element-deletion reduction
 * could achieve, since it doesn't charge anything for the rotation itself.
 *
 * For each pi in S_n (certified dist_n{n}.bin), for each of the n choices
 * of v, and (in the "with rotation" variant) each of the n-1 rotations of
 * the reduced array, look up d_{n-1} in dist_n{n-1}.bin and take
 * extra = d_n(pi) - d_{n-1}(...). Report:
 *   worst_no_rotation        = max_pi min_v extra   (natural array order only)
 *   worst_with_free_rotation = max_pi min_{v,k} extra (rotation cost ignored)
 * If worst_with_free_rotation > n-1 for some n, the whole "delete one
 * element" reduction family cannot give the target zero-slack recurrence,
 * even in the most generous (uncharged rotation) form.
 *
 * Usage: h3_reduction n tables_dir
 *   requires dist_n{n}.bin and dist_n{n-1}.bin (PROBLEM.md format) in
 *   tables_dir; n = 11, 12 tables are reproducible via exact/bfs_fast.c
 *   (not stored in git, see repo_state.md) but not required for n <= 10.
 * Version h3_reduction-1.0 (session 10, 14.09.2026).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 13
static uint64_t fact[MAXN + 1];
static void init_fact(int n) { fact[0] = 1; for (int i = 1; i <= n; i++) fact[i] = fact[i - 1] * (uint64_t)i; }

static uint64_t rank_perm(const uint8_t *p, int n) {
    uint64_t r = 0;
    for (int i = 0; i < n; i++) {
        int c = 0;
        for (int j = i + 1; j < n; j++) if (p[j] < p[i]) c++;
        r += (uint64_t)c * fact[n - 1 - i];
    }
    return r;
}

static uint8_t *load_table(const char *path, uint64_t total) {
    uint8_t *buf = malloc(total);
    FILE *f = fopen(path, "rb");
    if (!f || fread(buf, 1, total, f) != total) { fprintf(stderr, "cannot read %s\n", path); exit(1); }
    fclose(f);
    return buf;
}

static int next_perm(uint8_t *a, int n) {
    int i = n - 2;
    while (i >= 0 && a[i] >= a[i + 1]) i--;
    if (i < 0) return 0;
    int j = n - 1;
    while (a[j] <= a[i]) j--;
    uint8_t t = a[i]; a[i] = a[j]; a[j] = t;
    for (int l = i + 1, r = n - 1; l < r; l++, r--) { t = a[l]; a[l] = a[r]; a[r] = t; }
    return 1;
}

static void unrank_perm(uint64_t r, int n, uint8_t *out) {
    uint8_t avail[MAXN];
    for (int i = 0; i < n; i++) avail[i] = (uint8_t)i;
    int m = n;
    for (int i = n - 1; i >= 0; i--) {
        uint64_t f = fact[i];
        int idx = (int)(r / f);
        r %= f;
        out[n - 1 - i] = avail[idx];
        memmove(avail + idx, avail + idx + 1, (size_t)(m - idx - 1));
        m--;
    }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s n tables_dir\n", argv[0]); return 2; }
    int n = atoi(argv[1]);
    const char *dir = argv[2];
    int nm1 = n - 1;
    if (n < 5 || n > MAXN) { fprintf(stderr, "n out of range 5..%d\n", MAXN); return 2; }
    init_fact(n);
    uint64_t total_n = fact[n], total_nm1 = fact[nm1];

    char path[512];
    snprintf(path, sizeof path, "%s/dist_n%d.bin", dir, n);
    uint8_t *distn = load_table(path, total_n);
    snprintf(path, sizeof path, "%s/dist_n%d.bin", dir, nm1);
    uint8_t *distnm1 = load_table(path, total_nm1);

    uint8_t pi[MAXN], red[MAXN], rot[MAXN];
    for (int i = 0; i < n; i++) pi[i] = (uint8_t)i;

    int worst_norot = -1, worst_rot = -1;
    long argmax_norot = -1, argmax_rot = -1;
    long rank = -1;
    do {
        rank++;
        int dn = distn[rank_perm(pi, n)];
        int best_norot = 1000000, best_rot = 1000000;
        for (int v = 0; v < n; v++) {
            int m = 0;
            for (int i = 0; i < n; i++) {
                if (pi[i] == v) continue;
                red[m++] = pi[i] > v ? pi[i] - 1 : pi[i];
            }
            int dnm1 = distnm1[rank_perm(red, nm1)];
            int extra = dn - dnm1;
            if (extra < best_norot) best_norot = extra;
            if (extra < best_rot) best_rot = extra;
            memcpy(rot, red, nm1);
            for (int k = 1; k < nm1; k++) {
                uint8_t first = rot[0];
                for (int i = 0; i < nm1 - 1; i++) rot[i] = rot[i + 1];
                rot[nm1 - 1] = first;
                int dr = distnm1[rank_perm(rot, nm1)];
                int extra_r = dn - dr;
                if (extra_r < best_rot) best_rot = extra_r;
            }
        }
        if (best_norot > worst_norot) { worst_norot = best_norot; argmax_norot = rank; }
        if (best_rot > worst_rot) { worst_rot = best_rot; argmax_rot = rank; }
    } while (next_perm(pi, n));

    uint8_t out[MAXN];
    unrank_perm((uint64_t)argmax_rot, n, out);
    printf("n=%d target=%d worst_no_rotation=%d worst_with_free_rotation=%d rank_norot=%ld rank_rot=%ld argmax_rot=(",
           n, n - 1, worst_norot, worst_rot, argmax_norot, argmax_rot);
    for (int i = 0; i < n; i++) printf("%d%s", out[i], i + 1 < n ? "," : "");
    printf(") %s\n", worst_rot > n - 1 ? "EXCEEDS" : "ok");
    return 0;
}
