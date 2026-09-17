/* h13i_search.c — H13-I (докажи ли I(pi) <= floor((n-1)^2/4) конструктивно
 * через дешёвый кандидатный набор разрезов, а не полный перебор n^2 (q,c))?
 * Проверяет несколько кандидатов на роль "дешёвого сертификата" и считает,
 * сколько pi нарушают целевую границу, если ограничиться этим кандидатом:
 *
 *   local4  — худшее значение среди 2D-локальных минимумов inv(q,c)
 *             (4-окрестность (q +-1, c), (q, c +-1), тор)
 *   local8  — то же с 8-окрестностью (плюс диагонали)
 *   antipod — min по 4 фиксированным разрезам (0,0),(n/2,0),(0,n/2),(n/2,n/2)
 *             (только n чётно)
 *   qsweep  — min по q при фиксированном c = 0 (n кандидатов)
 *   L1/L2/Linf — c* минимизирует sum|t_i|, sum t_i^2, max|t_i| по
 *             t_i = centered((pi(i) - i + c) mod n), затем min по всем q
 *             при этом c* (n кандидатов)
 *
 * inv(w) считается перебором пар (n <= 11 — печально, но n <= 10 пробегается
 * за разумное время). I(pi) = min по всем n^2 разрезам — эталон для сравнения.
 *
 * Usage: h13i_search n
 * Version h13i_search-1.0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 11

static int inversions(const int *w, int n) {
    int inv = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (w[i] > w[j]) inv++;
    return inv;
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

static int inv_qc(const int *pi, int n, int q, int c, int *w) {
    int shift = q + 1 - c;
    for (int j = 0; j < n; j++) {
        int v = (pi[(q + 1 + j) % n] - shift) % n;
        if (v < 0) v += n;
        w[j] = v;
    }
    return inversions(w, n);
}

static int centered(int x, int n) { return (2 * x > n) ? x - n : x; }

int grid[MAXN][MAXN];

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s n\n", argv[0]); return 2; }
    int n = atoi(argv[1]);
    if (n < 3 || n > MAXN) { fprintf(stderr, "3 <= n <= %d\n", MAXN); return 2; }
    int bound = (n - 1) * (n - 1) / 4;
    int pi[MAXN], w[MAXN];
    for (int i = 0; i < n; i++) pi[i] = i;

    long total = 0;
    long viol_local4 = 0, viol_local8 = 0, viol_antipod = 0, viol_qsweep = 0;
    long viol_l1 = 0, viol_l2 = 0, viol_linf = 0;
    int worst_local4 = -1, worst_local8 = -1, worst_antipod = -1, worst_qsweep = -1;
    int worst_l1 = -1, worst_l2 = -1, worst_linf = -1;
    int global_max_I = -1;
    int half = n / 2;

    do {
        for (int q = 0; q < n; q++)
            for (int c = 0; c < n; c++)
                grid[q][c] = inv_qc(pi, n, q, c, w);

        int gmin = 999999;
        for (int q = 0; q < n; q++) for (int c = 0; c < n; c++) if (grid[q][c] < gmin) gmin = grid[q][c];
        if (gmin > global_max_I) global_max_I = gmin;

        /* local minima, 4- and 8-neighborhood (worst value among them) */
        int wl4 = -1, wl8 = -1;
        for (int q = 0; q < n; q++) for (int c = 0; c < n; c++) {
            int v = grid[q][c];
            int qm = (q - 1 + n) % n, qp = (q + 1) % n, cm = (c - 1 + n) % n, cp = (c + 1) % n;
            int is4 = (v <= grid[qm][c] && v <= grid[qp][c] && v <= grid[q][cm] && v <= grid[q][cp]);
            if (is4 && v > wl4) wl4 = v;
            int is8 = is4 && (v <= grid[qm][cm] && v <= grid[qm][cp] && v <= grid[qp][cm] && v <= grid[qp][cp]);
            if (is8 && v > wl8) wl8 = v;
        }
        if (wl4 > worst_local4) worst_local4 = wl4;
        if (wl4 > bound) viol_local4++;
        if (wl8 > worst_local8) worst_local8 = wl8;
        if (wl8 > bound) viol_local8++;

        if (n % 2 == 0) {
            int m = grid[0][0];
            if (grid[half][0] < m) m = grid[half][0];
            if (grid[0][half] < m) m = grid[0][half];
            if (grid[half][half] < m) m = grid[half][half];
            if (m > worst_antipod) worst_antipod = m;
            if (m > bound) viol_antipod++;
        }

        int mq = 999999;
        for (int q = 0; q < n; q++) if (grid[q][0] < mq) mq = grid[q][0];
        if (mq > worst_qsweep) worst_qsweep = mq;
        if (mq > bound) viol_qsweep++;

        int minsum1 = 999999, minsum2 = 999999, mininf = 999999, c1 = -1, c2 = -1, cinf = -1;
        for (int c = 0; c < n; c++) {
            int s1 = 0, s2 = 0, sinf = 0;
            for (int i = 0; i < n; i++) {
                int t = (pi[i] - i + c) % n; if (t < 0) t += n;
                int a = abs(centered(t, n));
                s1 += a; s2 += a * a; if (a > sinf) sinf = a;
            }
            if (s1 < minsum1) { minsum1 = s1; c1 = c; }
            if (s2 < minsum2) { minsum2 = s2; c2 = c; }
            if (sinf < mininf) { mininf = sinf; cinf = c; }
        }
        int b1 = 999999, b2 = 999999, binf = 999999;
        for (int q = 0; q < n; q++) if (grid[q][c1] < b1) b1 = grid[q][c1];
        for (int q = 0; q < n; q++) if (grid[q][c2] < b2) b2 = grid[q][c2];
        for (int q = 0; q < n; q++) if (grid[q][cinf] < binf) binf = grid[q][cinf];
        if (b1 > worst_l1) worst_l1 = b1;
        if (b1 > bound) viol_l1++;
        if (b2 > worst_l2) worst_l2 = b2;
        if (b2 > bound) viol_l2++;
        if (binf > worst_linf) worst_linf = binf;
        if (binf > bound) viol_linf++;

        total++;
    } while (next_perm(pi, n));

    printf("{\"n\": %d, \"bound\": %d, \"total\": %ld, \"global_max_I\": %d,\n", n, bound, total, global_max_I);
    printf(" \"local4\": {\"worst\": %d, \"violations\": %ld},\n", worst_local4, viol_local4);
    printf(" \"local8\": {\"worst\": %d, \"violations\": %ld},\n", worst_local8, viol_local8);
    printf(" \"antipod4pt\": {\"worst\": %d, \"violations\": %ld},\n", worst_antipod, viol_antipod);
    printf(" \"qsweep_c0\": {\"worst\": %d, \"violations\": %ld},\n", worst_qsweep, viol_qsweep);
    printf(" \"L1_center_c\": {\"worst\": %d, \"violations\": %ld},\n", worst_l1, viol_l1);
    printf(" \"L2_center_c\": {\"worst\": %d, \"violations\": %ld},\n", worst_l2, viol_l2);
    printf(" \"Linf_center_c\": {\"worst\": %d, \"violations\": %ld}}\n", worst_linf, viol_linf);
    return 0;
}
