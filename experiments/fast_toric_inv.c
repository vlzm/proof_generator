/* fast_toric_inv.c — O(n) per rotation ( O(n^2) per permutation) algorithm for
 * I(pi) = min over double cuts (q, c) of inv(w), the same quantity as
 * experiments/line_profile.c (which computes it in O(n^4) per permutation by
 * brute force). Used to test H13-I (docs/notes/h13_line_model.md §6):
 * I(pi) <= floor((n-1)^2/4) for all pi.
 *
 * Recurrence used (see docs/notes/h13_line_model.md session 9 write-up):
 * for a fixed rotation v of pi (position cut), let pos = inverse permutation
 * of v. As the value-cut shifts s -> s+1, inversions change by exactly
 * n-1-2*pos[s] (the element with original value s jumps from rank 0 to rank
 * n-1). So min over s of inversions = inv(v) + min_{t=0..n-1} prefix_sum(a),
 * a_k = n-1-2*pos[k], prefix_sum(0) = 0. And rotating positions (v -> v shifted
 * left by one, dropping v[0] to the back) changes inv(v) by exactly
 * n-1-2*v[0] (v[0] is a value in 0..n-1, so the count of smaller elements
 * among the rest is exactly v[0]). Both updates are O(n) (pos array shifts by
 * -1 mod n on rotation), giving O(n) per q and O(n^2) per permutation total,
 * versus O(n^4) for the direct approach.
 *
 * C_q(pi) := min_c inv(w_{q,c}) (fast_I minimised over a single q instead of
 * all q) satisfies C_q = inv(v) + min_t S_t as above. H13-I-avg (session 9,
 * `docs/notes/h13_line_model.md` §7): (1/n) sum_q C_q(pi) <= floor((n-1)^2/4)
 * for every pi (would imply H13-I, since min_q C_q <= average); verified
 * exhaustively 4 <= n <= 12, tight exactly on reflections. Modes q0exhaust
 * and pairexhaust test (and refute, for n large enough) two candidate
 * universal single/paired choices of q that were checked and rejected as
 * insufficient (do not repeat).
 *
 * Usage:
 *   fast_toric_inv validate n   — compare fast and brute-force I(pi) for all
 *     permutations of size n (n <= 9), abort with the first mismatch.
 *   fast_toric_inv exhaust n    — exhaustive over all n! permutations
 *     (Lehmer/lex order), report max I(pi), argmax, count at each I value,
 *     and check against floor((n-1)^2/4).
 *   fast_toric_inv sample n trials seed [family] — random or structured
 *     samples (family: random | affine | reflection | rotsigma), report max I.
 *   fast_toric_inv localsearch n restarts stall seed — hill-climb (strict
 *     improvement, single transpositions) trying to push I(pi) above the
 *     bound; adversarial counterexample search.
 *   fast_toric_inv avgexhaust n — exhaustive: max_pi (1/n) sum_q C_q(pi),
 *     the H13-I-avg lemma above.
 *   fast_toric_inv q0exhaust n — exhaustive: max_pi C_0(pi) (q fixed at 0
 *     only); known to exceed the bound starting at n = 7 (REFUTED as a
 *     universal single-q rule).
 *   fast_toric_inv pairexhaust n off — exhaustive: max_pi min(C_0, C_off)
 *     (a single fixed pairing offset); known to fail for some offsets and
 *     n (e.g. off=1 at n=10,11; off=n/2 at n=8) — REFUTED as a universal
 *     rule, the good q genuinely depends on pi.
 * Version fast_toric_inv-1.1.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 4096

static int brute_I(const int *pi, int n) {
    int best = 1 << 30;
    int w[MAXN];
    for (int q = 0; q < n; q++) {
        for (int c = 0; c < n; c++) {
            int shift = ((q + 1 - c) % n + n) % n;
            for (int j = 0; j < n; j++) {
                int v = (pi[(q + 1 + j) % n] - shift) % n;
                if (v < 0) v += n;
                w[j] = v;
            }
            int inv = 0;
            for (int i = 0; i < n; i++)
                for (int j2 = i + 1; j2 < n; j2++)
                    if (w[i] > w[j2]) inv++;
            if (inv < best) best = inv;
        }
    }
    return best;
}

/* Fast I(pi): O(n^2). */
static int fast_I(const int *pi, int n) {
    int v[MAXN], pos[MAXN];
    for (int j = 0; j < n; j++) v[j] = pi[j];
    /* initial inversions of v, O(n^2) (fine; dominated by outer n loop anyway) */
    int inv0 = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (v[i] > v[j]) inv0++;
    for (int j = 0; j < n; j++) pos[v[j]] = j;

    int best = 1 << 30;
    for (int q = 0; q < n; q++) {
        /* min over s of inv0 + prefix_sum(a), a_k = n-1-2*pos[k] */
        long acc = 0, minacc = 0;
        for (int k = 0; k < n; k++) {
            long a = (long)(n - 1) - 2L * pos[k];
            acc += a;
            if (acc < minacc) minacc = acc;
        }
        long cand = (long)inv0 + minacc;
        if (cand < best) best = (int)cand;

        if (q == n - 1) break;
        /* rotate v left by one: drop v[0], append at end */
        int dropped = v[0];
        inv0 += (n - 1 - 2 * dropped);
        for (int j = 0; j < n - 1; j++) v[j] = v[j + 1];
        v[n - 1] = dropped;
        /* update pos: every value's position shifts by -1 mod n */
        for (int k = 0; k < n; k++) pos[k] = (pos[k] - 1 + n) % n;
    }
    return best;
}

/* Fill Cq[0..n-1] with C_q(pi) = min_c inv(w_{q,c}) for every q. */
static void all_Cq(const int *pi, int n, int *Cq) {
    int v[MAXN], pos[MAXN];
    for (int j = 0; j < n; j++) v[j] = pi[j];
    int inv0 = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (v[i] > v[j]) inv0++;
    for (int j = 0; j < n; j++) pos[v[j]] = j;
    for (int q = 0; q < n; q++) {
        long acc = 0, minacc = 0;
        for (int k = 0; k < n; k++) {
            long a = (long)(n - 1) - 2L * pos[k];
            acc += a;
            if (acc < minacc) minacc = acc;
        }
        Cq[q] = inv0 + (int)minacc;
        if (q == n - 1) break;
        int dropped = v[0];
        inv0 += (n - 1 - 2 * dropped);
        for (int j = 0; j < n - 1; j++) v[j] = v[j + 1];
        v[n - 1] = dropped;
        for (int k = 0; k < n; k++) pos[k] = (pos[k] - 1 + n) % n;
    }
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

static void run_validate(int n) {
    if (n > 9) { fprintf(stderr, "validate: n <= 9 only (brute force too slow)\n"); exit(2); }
    int pi[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;
    long rank = 0, mism = 0;
    do {
        int a = fast_I(pi, n);
        int b = brute_I(pi, n);
        if (a != b) {
            printf("MISMATCH rank=%ld fast=%d brute=%d pi=(", rank, a, b);
            for (int i = 0; i < n; i++) printf("%d%s", pi[i], i + 1 < n ? "," : "");
            printf(")\n");
            mism++;
        }
        rank++;
    } while (next_perm(pi, n));
    printf("{\"mode\":\"validate\",\"n\":%d,\"count\":%ld,\"mismatches\":%ld}\n", n, rank, mism);
}

static void run_exhaust(int n) {
    if (n > 13) { fprintf(stderr, "exhaust: n <= 13 (time)\n"); exit(2); }
    int pi[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;
    long total = 1;
    for (int k = 2; k <= n; k++) total *= k;
    int bound = ((n - 1) * (n - 1)) / 4;
    long cnt_at[4096]; memset(cnt_at, 0, sizeof cnt_at);
    int maxI = -1; long argmax = -1, nviol = 0; long firstviol = -1;
    long rank = 0;
    do {
        int I = fast_I(pi, n);
        cnt_at[I]++;
        if (I > maxI) { maxI = I; argmax = rank; }
        if (I > bound) { nviol++; if (firstviol < 0) firstviol = rank; }
        rank++;
    } while (next_perm(pi, n));
    printf("{\"mode\":\"exhaust\",\"n\":%d,\"count\":%ld,\"bound_floor((n-1)^2/4)\":%d,"
           "\"max_I\":%d,\"argmax_rank\":%ld,\"violations\":%ld,\"first_violation_rank\":%ld,\"by_I\":[",
           n, total, bound, maxI, argmax, nviol, firstviol);
    int first = 1;
    for (int i = 0; i < 4096; i++) if (cnt_at[i]) {
        printf("%s[%d,%ld]", first ? "" : ",", i, cnt_at[i]);
        first = 0;
    }
    printf("]}\n");
}

static int gcd_local(int a, int b) { while (b) { int t = b; b = a % b; a = t; } return a; }

static unsigned long rngstate;
static unsigned long xrand(void) {
    rngstate ^= rngstate << 13; rngstate ^= rngstate >> 7; rngstate ^= rngstate << 17;
    return rngstate;
}

static void run_sample(int n, long trials, unsigned long seed, const char *family) {
    if (n > MAXN) { fprintf(stderr, "n too large\n"); exit(2); }
    rngstate = seed ? seed : 88172645463325252UL;
    int pi[MAXN];
    int bound = ((n - 1) * (n - 1)) / 4;
    int maxI = -1; long argbest = -1; long nviol = 0;
    long t;
    for (t = 0; t < trials; t++) {
        if (strcmp(family, "random") == 0) {
            for (int i = 0; i < n; i++) pi[i] = i;
            for (int i = n - 1; i > 0; i--) {
                int j = xrand() % (i + 1);
                int tmp = pi[i]; pi[i] = pi[j]; pi[j] = tmp;
            }
        } else if (strcmp(family, "affine") == 0) {
            int a, b;
            do { a = 1 + (int)(xrand() % (n - 1)); } while (gcd_local(a, n) != 1);
            b = (int)(xrand() % n);
            for (int i = 0; i < n; i++) pi[i] = ((long)a * i + b) % n;
        } else if (strcmp(family, "reflection") == 0) {
            int h = (int)(xrand() % n);
            for (int i = 0; i < n; i++) { int v = h - i; v %= n; if (v < 0) v += n; pi[i] = v; }
        } else if (strcmp(family, "rotsigma") == 0) {
            int k = (int)(xrand() % n);
            for (int i = 0; i < n; i++) {
                int base = n - 1 - i; /* sigma_n(i) = n-1-i */
                pi[i] = ((base + k) % n + n) % n;
            }
        } else {
            fprintf(stderr, "unknown family %s\n", family);
            exit(2);
        }
        int I = fast_I(pi, n);
        if (I > maxI) { maxI = I; argbest = t; }
        if (I > bound) nviol++;
    }
    printf("{\"mode\":\"sample\",\"family\":\"%s\",\"n\":%d,\"trials\":%ld,\"seed\":%lu,"
           "\"bound_floor((n-1)^2/4)\":%d,\"max_I\":%d,\"argbest_trial\":%ld,\"violations\":%ld}\n",
           family, n, trials, seed, bound, maxI, argbest, nviol);
}

/* Local search: hill-climb pi under random transpositions to try to push
 * I(pi) above floor((n-1)^2/4). Multiple restarts, each restart runs until
 * `stall` consecutive failed proposals. Reports the best I found overall and
 * whether it beats the bound (adversarial search for a counterexample,
 * complementing the structured families above). */
static void run_localsearch(int n, long restarts, long stall, unsigned long seed) {
    rngstate = seed ? seed : 88172645463325252UL;
    int pi[MAXN], best_pi[MAXN];
    int bound = ((n - 1) * (n - 1)) / 4;
    int globalbest = -1;
    for (long r = 0; r < restarts; r++) {
        for (int i = 0; i < n; i++) pi[i] = i;
        for (int i = n - 1; i > 0; i--) {
            int j = xrand() % (i + 1);
            int tmp = pi[i]; pi[i] = pi[j]; pi[j] = tmp;
        }
        int cur = fast_I(pi, n);
        long fails = 0;
        while (fails < stall) {
            int i = xrand() % n, j = xrand() % n;
            if (i == j) continue;
            int tmp = pi[i]; pi[i] = pi[j]; pi[j] = tmp;
            int cand = fast_I(pi, n);
            if (cand > cur) { cur = cand; fails = 0; }
            else { pi[j] = pi[i]; pi[i] = tmp; fails++; }
        }
        if (cur > globalbest) {
            globalbest = cur;
            for (int i = 0; i < n; i++) best_pi[i] = pi[i];
        }
    }
    printf("{\"mode\":\"localsearch\",\"n\":%d,\"restarts\":%ld,\"stall\":%ld,\"seed\":%lu,"
           "\"bound_floor((n-1)^2/4)\":%d,\"max_I\":%d,\"violation\":%s,\"best_pi\":[",
           n, restarts, stall, seed, bound, globalbest, globalbest > bound ? "true" : "false");
    for (int i = 0; i < n; i++) printf("%s%d", i ? "," : "", best_pi[i]);
    printf("]}\n");
}

static void run_avgexhaust(int n) {
    if (n > 13) { fprintf(stderr, "avgexhaust: n <= 13 (time)\n"); exit(2); }
    int pi[MAXN], Cq[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;
    long total = 1;
    for (int k = 2; k <= n; k++) total *= k;
    long bound_n = (long)(((n - 1) * (n - 1)) / 4) * n;
    long cnt = 0, viol = 0, maxsum = -1, argmax = -1;
    do {
        all_Cq(pi, n, Cq);
        long s = 0;
        for (int q = 0; q < n; q++) s += Cq[q];
        if (s > maxsum) { maxsum = s; argmax = cnt; }
        if (s > bound_n) viol++;
        cnt++;
    } while (next_perm(pi, n));
    printf("{\"mode\":\"avgexhaust\",\"n\":%d,\"count\":%ld,\"n_times_bound\":%ld,"
           "\"max_sum_Cq\":%ld,\"max_avg\":%.6f,\"bound\":%d,\"argmax_rank\":%ld,\"violations\":%ld}\n",
           n, cnt, bound_n, maxsum, (double)maxsum / n, (n - 1) * (n - 1) / 4, argmax, viol);
}

static void run_q0exhaust(int n) {
    if (n > 13) { fprintf(stderr, "q0exhaust: n <= 13 (time)\n"); exit(2); }
    int pi[MAXN], Cq[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;
    long cnt = 0, viol = 0; int bound = (n - 1) * (n - 1) / 4; int maxv = -1;
    do {
        all_Cq(pi, n, Cq);
        if (Cq[0] > maxv) maxv = Cq[0];
        if (Cq[0] > bound) viol++;
        cnt++;
    } while (next_perm(pi, n));
    printf("{\"mode\":\"q0exhaust\",\"n\":%d,\"count\":%ld,\"bound\":%d,\"max_C0\":%d,\"violations\":%ld}\n",
           n, cnt, bound, maxv, viol);
}

static void run_pairexhaust(int n, int off) {
    if (n > 13) { fprintf(stderr, "pairexhaust: n <= 13 (time)\n"); exit(2); }
    int pi[MAXN], Cq[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;
    long cnt = 0, viol = 0; int bound = (n - 1) * (n - 1) / 4; int maxv = -1;
    int j = ((off % n) + n) % n;
    do {
        all_Cq(pi, n, Cq);
        int m = Cq[0] < Cq[j] ? Cq[0] : Cq[j];
        if (m > maxv) maxv = m;
        if (m > bound) viol++;
        cnt++;
    } while (next_perm(pi, n));
    printf("{\"mode\":\"pairexhaust\",\"n\":%d,\"off\":%d,\"count\":%ld,\"bound\":%d,"
           "\"max_min_C0_Coff\":%d,\"violations\":%ld}\n",
           n, off, cnt, bound, maxv, viol);
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s validate n | exhaust n | sample n trials seed [family]\n", argv[0]);
        return 2;
    }
    const char *mode = argv[1];
    int n = atoi(argv[2]);
    if (strcmp(mode, "validate") == 0) run_validate(n);
    else if (strcmp(mode, "exhaust") == 0) run_exhaust(n);
    else if (strcmp(mode, "sample") == 0) {
        if (argc < 5) { fprintf(stderr, "sample needs trials seed [family]\n"); return 2; }
        long trials = atol(argv[3]);
        unsigned long seed = strtoul(argv[4], NULL, 10);
        const char *family = argc >= 6 ? argv[5] : "random";
        run_sample(n, trials, seed, family);
    } else if (strcmp(mode, "localsearch") == 0) {
        if (argc < 6) { fprintf(stderr, "localsearch needs restarts stall seed\n"); return 2; }
        long restarts = atol(argv[3]);
        long stall = atol(argv[4]);
        unsigned long seed = strtoul(argv[5], NULL, 10);
        run_localsearch(n, restarts, stall, seed);
    } else if (strcmp(mode, "avgexhaust") == 0) run_avgexhaust(n);
    else if (strcmp(mode, "q0exhaust") == 0) run_q0exhaust(n);
    else if (strcmp(mode, "pairexhaust") == 0) {
        if (argc < 4) { fprintf(stderr, "pairexhaust needs off\n"); return 2; }
        run_pairexhaust(n, atoi(argv[3]));
    } else {
        fprintf(stderr, "unknown mode %s\n", mode);
        return 2;
    }
    return 0;
}
