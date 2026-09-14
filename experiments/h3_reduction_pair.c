/* h3_reduction_pair.c -- two-level variant of h3_reduction.c: test whether
 * removing a PAIR of array elements at once (v1, v2 out of n choose 2,
 * array closes up, remaining n-2 values relabelled to 0..n-3) can satisfy
 * the two-step telescoped H3 budget d_n(pi) - d_{n-2}(reduced) <= 2n-3
 * (= (n-1)+(n-2), matching two levels of the H3 recurrence T_n <=
 * T_{n-1}+(n-1)) for every pi in S_n. This is strictly more generous than
 * composing two single-element removals (h3_reduction.c), since it
 * searches jointly over both elements removed together rather than one
 * at a time -- so if it also fails, no k=1 or k=2 "delete elements,
 * relabel" scheme works, at least not with a fixed small k.
 *
 * No rotation freedom here (h3_reduction.c already showed rotation only
 * delays failure, not prevents it, and pairs already multiply the search
 * space by ~n/2 over singles).
 *
 * Usage: h3_reduction_pair n tables_dir
 *   requires dist_n{n}.bin and dist_n{n-2}.bin (PROBLEM.md format).
 * Version h3_reduction_pair-1.0 (session 11, 14.09.2026).
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
    int nm2 = n - 2;
    if (n < 6 || n > MAXN) { fprintf(stderr, "n out of range 6..%d\n", MAXN); return 2; }
    init_fact(n);
    uint64_t total_n = fact[n], total_nm2 = fact[nm2];

    char path[512];
    snprintf(path, sizeof path, "%s/dist_n%d.bin", dir, n);
    uint8_t *distn = load_table(path, total_n);
    snprintf(path, sizeof path, "%s/dist_n%d.bin", dir, nm2);
    uint8_t *distnm2 = load_table(path, total_nm2);

    uint8_t pi[MAXN], red[MAXN];
    for (int i = 0; i < n; i++) pi[i] = (uint8_t)i;

    int target = 2 * n - 3;
    int worst = -1;
    long argmax = -1;
    long rank = -1;
    do {
        rank++;
        int dn = distn[rank_perm(pi, n)];
        int best = 1000000;
        for (int v1 = 0; v1 < n; v1++) {
            for (int v2 = v1 + 1; v2 < n; v2++) {
                int m = 0;
                for (int i = 0; i < n; i++) {
                    int x = pi[i];
                    if (x == v1 || x == v2) continue;
                    int shift = (x > v1) + (x > v2);
                    red[m++] = (uint8_t)(x - shift);
                }
                int dnm2 = distnm2[rank_perm(red, nm2)];
                int extra = dn - dnm2;
                if (extra < best) best = extra;
            }
        }
        if (best > worst) { worst = best; argmax = rank; }
    } while (next_perm(pi, n));

    uint8_t out[MAXN];
    unrank_perm((uint64_t)argmax, n, out);
    printf("n=%d target=2n-3=%d worst=%d rank=%ld argmax=(", n, target, worst, argmax);
    for (int i = 0; i < n; i++) printf("%d%s", out[i], i + 1 < n ? "," : "");
    printf(") %s\n", worst > target ? "EXCEEDS" : "ok");
    return 0;
}
