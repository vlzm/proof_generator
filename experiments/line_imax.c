/* line_imax.c -- exhaustive test of H13-I:  I(pi) <= floor((n-1)^2/4)
 * for every toroidal permutation class of Z_n.
 *
 * I(pi) = min over the n^2 double cuts (a, b) of inv(w), w_j = (pi(a+j) - b) mod n
 * (the convention of experiments/line_model.py with a = q+1, b = q+1-c).
 *
 * I(pi) is invariant under pi(i) -> pi(i + s) - t, so it is enough to run over the
 * (n-1)! permutations with pi(0) = 0: every toroidal class has such a representative
 * (in fact exactly n of them).  For one pi the whole n x n table of inv is built in
 * O(n^2) from the two increment rules
 *     inv(a+1, b) - inv(a, b) = (n-1) - 2 * ((pi(a) - b) mod n),
 *     inv(a, b+1) - inv(a, b) = (n-1) - 2 * ((posval(b) - a) mod n),
 * where posval(b) is the position carrying the value b; only inv(0,0) is computed
 * by the O(n^2) double loop.
 *
 * Output: max I, the number of classes attaining it, the first attaining pi, and the
 * number of violations of floor((n-1)^2/4).  Progress lines every CHECKPOINT
 * permutations so a long run can be monitored and stopped.
 *
 * Usage: line_imax n [checkpoint]
 * Version line_imax-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAXN 14

static int n;
static int pi_[MAXN], posval[MAXN];
static int row[MAXN];          /* inv(a, b) for the current a */

/* inversions of the line at cut (0,0) */
static int inv00(void) {
    int c = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (pi_[i] > pi_[j]) c++;
    return c;
}

static int I_of_pi(void) {
    for (int i = 0; i < n; i++) posval[pi_[i]] = i;
    int best = 1 << 30;
    int cur = inv00();                      /* inv(0,0) */
    for (int a = 0; a < n; a++) {
        /* fill the row b = 0..n-1 starting from inv(a,0) = cur */
        int v = cur;
        for (int b = 0; b < n; b++) {
            if (v < best) best = v;
            /* step b -> b+1 */
            int p = posval[b] - a; if (p < 0) p += n;
            v += (n - 1) - 2 * p;
        }
        /* step a -> a+1 at b = 0 */
        int r = pi_[a];                     /* (pi(a) - 0) mod n */
        cur += (n - 1) - 2 * r;
    }
    return best;
}

static int next_perm_tail(void) {           /* next permutation of pi_[1..n-1] */
    int i = n - 2;
    while (i >= 1 && pi_[i] >= pi_[i + 1]) i--;
    if (i < 1) return 0;
    int j = n - 1;
    while (pi_[j] <= pi_[i]) j--;
    int t = pi_[i]; pi_[i] = pi_[j]; pi_[j] = t;
    for (int l = i + 1, r = n - 1; l < r; l++, r--) { t = pi_[l]; pi_[l] = pi_[r]; pi_[r] = t; }
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s n [checkpoint]\n", argv[0]); return 2; }
    n = atoi(argv[1]);
    if (n < 3 || n > MAXN) { fprintf(stderr, "n out of range\n"); return 2; }
    long long checkpoint = (argc >= 3) ? atoll(argv[2]) : 100000000LL;
    int bound = ((n - 1) * (n - 1)) / 4;

    for (int i = 0; i < n; i++) pi_[i] = i;
    long long cnt = 0, viol = 0, attain = 0;
    int maxI = -1;
    int first[MAXN];
    clock_t t0 = clock();
    do {
        int I = I_of_pi();
        if (I > maxI) { maxI = I; attain = 0; }
        if (I == maxI) { if (attain == 0) for (int i = 0; i < n; i++) first[i] = pi_[i]; attain++; }
        if (I > bound) {
            viol++;
            if (viol <= 5) {
                printf("VIOLATION I=%d > %d  pi =", I, bound);
                for (int i = 0; i < n; i++) printf(" %d", pi_[i]);
                printf("\n"); fflush(stdout);
            }
        }
        cnt++;
        if (cnt % checkpoint == 0) {
            printf("# progress %lld  maxI=%d  viol=%lld  %.1fs\n",
                   cnt, maxI, viol, (double)(clock() - t0) / CLOCKS_PER_SEC);
            fflush(stdout);
        }
    } while (next_perm_tail());

    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("{\"version\": \"line_imax-1.0\", \"n\": %d, \"representatives\": %lld, "
           "\"bound_floor_(n-1)^2/4\": %d, \"max_I\": %d, \"argmax_count\": %lld, "
           "\"violations\": %lld, \"seconds\": %.2f, \"first_argmax\": [", n, cnt, bound, maxI, attain, viol, secs);
    for (int i = 0; i < n; i++) printf("%s%d", i ? ", " : "", first[i]);
    printf("]}\n");
    return viol ? 1 : 0;
}
