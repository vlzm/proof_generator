/* Fast BFS + certificate for the LRX Cayley graph.
 *
 * Same specification as exact/bfs.py and oracle/verify_distances.py:
 * - permutations of 0..n-1, Lehmer rank indexing (count-smaller-to-right);
 * - moves: L = left cyclic shift, R = right cyclic shift, X = swap first two;
 * - table: one byte per state, 0xFF = unvisited.
 *
 * Must agree byte-for-byte with the Python implementation on small n
 * (cross-checked externally); the certificate re-checks the table
 * independently of the BFS order.
 *
 * Usage: bfs_fast n out.bin        — run BFS from id_n, write table
 *        bfs_fast n table.bin cert — certify an existing table
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

static uint64_t fact[21];

static void init_fact(int n) {
    fact[0] = 1;
    for (int i = 1; i <= n; i++) fact[i] = fact[i - 1] * (uint64_t)i;
}

static uint64_t rank_perm(const uint8_t *p, int n) {
    uint64_t r = 0;
    for (int i = 0; i < n; i++) {
        int c = 0;
        for (int j = i + 1; j < n; j++) if (p[j] < p[i]) c++;
        r += (uint64_t)c * fact[n - 1 - i];
    }
    return r;
}

static void unrank_perm(uint64_t r, int n, uint8_t *out) {
    uint8_t avail[21];
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

static int run_bfs(int n, const char *out_path) {
    init_fact(n);
    uint64_t total = fact[n];
    uint8_t *dist = malloc(total);
    if (!dist) { fprintf(stderr, "alloc dist failed\n"); return 1; }
    memset(dist, 0xFF, total);

    uint32_t *frontier = malloc(sizeof(uint32_t) * (total / 4 + 16));
    uint32_t *next = malloc(sizeof(uint32_t) * (total / 4 + 16));
    if (!frontier || !next) { fprintf(stderr, "alloc frontier failed\n"); return 1; }

    uint8_t p[21], q[21];
    for (int i = 0; i < n; i++) p[i] = (uint8_t)i;
    uint64_t r0 = rank_perm(p, n);
    dist[r0] = 0;
    frontier[0] = (uint32_t)r0;
    uint64_t flen = 1, visited = 1;
    int d = 0;
    time_t t0 = time(NULL);

    while (flen > 0) {
        d++;
        if (d >= 0xFF) { fprintf(stderr, "distance overflow\n"); return 1; }
        uint64_t nlen = 0;
        for (uint64_t k = 0; k < flen; k++) {
            uint64_t r = frontier[k];
            unrank_perm(r, n, p);
            /* L */
            memcpy(q, p + 1, (size_t)(n - 1)); q[n - 1] = p[0];
            uint64_t rq = rank_perm(q, n);
            if (dist[rq] == 0xFF) { dist[rq] = (uint8_t)d; next[nlen++] = (uint32_t)rq; }
            /* R */
            q[0] = p[n - 1]; memcpy(q + 1, p, (size_t)(n - 1));
            rq = rank_perm(q, n);
            if (dist[rq] == 0xFF) { dist[rq] = (uint8_t)d; next[nlen++] = (uint32_t)rq; }
            /* X */
            memcpy(q, p, (size_t)n); q[0] = p[1]; q[1] = p[0];
            rq = rank_perm(q, n);
            if (dist[rq] == 0xFF) { dist[rq] = (uint8_t)d; next[nlen++] = (uint32_t)rq; }
        }
        visited += nlen;
        fprintf(stderr, "layer %d: %llu new, %llu/%llu visited, %lds\n",
                d, (unsigned long long)nlen, (unsigned long long)visited,
                (unsigned long long)total, (long)(time(NULL) - t0));
        uint32_t *tmp = frontier; frontier = next; next = tmp;
        flen = nlen;
    }
    int maxd = 0;
    for (uint64_t r = 0; r < total; r++) {
        if (dist[r] == 0xFF) { fprintf(stderr, "INCOMPLETE at %llu\n",
                                       (unsigned long long)r); return 1; }
        if (dist[r] > maxd) maxd = dist[r];
    }
    FILE *f = fopen(out_path, "wb");
    if (!f || fwrite(dist, 1, total, f) != total) {
        fprintf(stderr, "write failed\n"); return 1;
    }
    fclose(f);
    printf("n=%d complete, D_n=%d, states=%llu\n", n, maxd,
           (unsigned long long)total);
    return 0;
}

static int certify(int n, const char *path) {
    init_fact(n);
    uint64_t total = fact[n];
    uint8_t *dist = malloc(total);
    FILE *f = fopen(path, "rb");
    if (!f || fread(dist, 1, total, f) != total) {
        fprintf(stderr, "read failed\n"); return 1;
    }
    fclose(f);

    uint8_t p[21], q[21];
    for (int i = 0; i < n; i++) p[i] = (uint8_t)i;
    uint64_t id_rank = rank_perm(p, n);
    uint64_t zeros = 0;
    for (uint64_t r = 0; r < total; r++) {
        if (dist[r] == 0xFF) { printf("FAIL: unvisited %llu\n",
                                      (unsigned long long)r); return 1; }
        if (dist[r] == 0) {
            zeros++;
            if (r != id_rank) { printf("FAIL: zero not at id\n"); return 1; }
        }
    }
    if (zeros != 1) { printf("FAIL: %llu zeros\n", (unsigned long long)zeros); return 1; }

    time_t t0 = time(NULL);
    for (uint64_t r = 0; r < total; r++) {
        int h = dist[r];
        unrank_perm(r, n, p);
        int has_down = 0;
        for (int m = 0; m < 3; m++) {
            if (m == 0) { memcpy(q, p + 1, (size_t)(n - 1)); q[n - 1] = p[0]; }
            else if (m == 1) { q[0] = p[n - 1]; memcpy(q + 1, p, (size_t)(n - 1)); }
            else { memcpy(q, p, (size_t)n); q[0] = p[1]; q[1] = p[0]; }
            int hq = dist[rank_perm(q, n)];
            int diff = h - hq;
            if (diff < -1 || diff > 1) {
                printf("FAIL: edge violation at %llu (%d vs %d)\n",
                       (unsigned long long)r, h, hq);
                return 1;
            }
            if (hq == h - 1) has_down = 1;
        }
        if (h > 0 && !has_down) {
            printf("FAIL: no descending neighbour at %llu\n",
                   (unsigned long long)r);
            return 1;
        }
        if ((r & 0xFFFFFFF) == 0 && r)
            fprintf(stderr, "certified %llu/%llu, %lds\n",
                    (unsigned long long)r, (unsigned long long)total,
                    (long)(time(NULL) - t0));
    }
    printf("n=%d certificate PASS (all %llu states, all edges)\n", n,
           (unsigned long long)total);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s n out.bin [cert]\n", argv[0]); return 1; }
    int n = atoi(argv[1]);
    if (n < 2 || n > 13) { fprintf(stderr, "n out of range\n"); return 1; }
    if (argc >= 4 && strcmp(argv[3], "cert") == 0) return certify(n, argv[2]);
    return run_bfs(n, argv[2]);
}
