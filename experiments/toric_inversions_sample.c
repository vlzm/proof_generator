/*
 * Structured + random sampling of I(pi) for n beyond the exhaustive range
 * of toric_inversions_scan.c (H13-I / C33). Uses the same O(n^2) gradient
 * identity (C37) per permutation, so n up to a few hundred is cheap.
 *
 * Families sampled per n: all n reflections pi(i) = (h - i) mod n; affine
 * maps pi(i) = (a*i + b) mod n for a coprime with n; uniform random
 * permutations; reflections perturbed by 1-3 random transpositions.
 *
 * Goal: look for any pi with I(pi) > floor((n-1)^2/4) at n where exhaustive
 * search is infeasible (n >= 13), per AGENTS.md rule 9 (sampling ladder for
 * large n). This is a search for counterexamples, not a proof.
 *
 * Usage: ./toric_inversions_sample [n1 n2 ...]  (defaults to a built-in list)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int inv_count(const int *a, int n) {
    long long c = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (a[i] > a[j]) c++;
    return (int)c;
}

/* O(n^2) via the C37 gradient identity; O(n) extra memory. */
static long long compute_I(const int *pi, int n, const int *invpi) {
    long long *row = malloc(sizeof(long long) * n);
    long long g00 = inv_count(pi, n);
    long long mn = g00;
    long long ga0 = g00;
    for (int a = 0; a < n; a++) {
        row[0] = ga0;
        if (row[0] < mn) mn = row[0];
        for (int b = 0; b < n - 1; b++) {
            int jstar = invpi[b] - a;
            jstar %= n; if (jstar < 0) jstar += n;
            row[b + 1] = row[b] + (n - 1 - 2 * jstar);
            if (row[b + 1] < mn) mn = row[b + 1];
        }
        if (a < n - 1) ga0 += n - 1 - 2 * pi[a];
    }
    free(row);
    return mn;
}

static void fisher_yates(int *a, int n) {
    for (int i = n - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        int t = a[i]; a[i] = a[j]; a[j] = t;
    }
}

static void scan_n(int n, int random_trials) {
    int *pi_ = malloc(sizeof(int) * n);
    int *invpi_ = malloc(sizeof(int) * n);
    long long bound = (long long)(n - 1) * (n - 1) / 4;
    long long worst = -1;
    int worst_kind = -1; /* 0=reflection 1=affine 2=random 3=perturbed reflection */

    for (int h = 0; h < n; h++) {
        for (int i = 0; i < n; i++) { int v = ((h - i) % n + n) % n; pi_[i] = v; }
        for (int i = 0; i < n; i++) invpi_[pi_[i]] = i;
        long long val = compute_I(pi_, n, invpi_);
        if (val > worst) { worst = val; worst_kind = 0; }
    }

    for (int a = 2; a < n && a < 60; a++) {
        int step_b = n > 40 ? (n / 20) + 1 : 1;
        for (int b = 0; b < n; b += step_b) {
            int *chk = calloc(n, sizeof(int));
            int ok = 1;
            for (int i = 0; i < n; i++) {
                int v = ((long long)a * i + b) % n; if (v < 0) v += n;
                pi_[i] = v;
                if (chk[v]) { ok = 0; }
                chk[v] = 1;
            }
            free(chk);
            if (!ok) continue;
            for (int i = 0; i < n; i++) invpi_[pi_[i]] = i;
            long long val = compute_I(pi_, n, invpi_);
            if (val > worst) { worst = val; worst_kind = 1; }
        }
    }

    for (int t = 0; t < random_trials; t++) {
        for (int i = 0; i < n; i++) pi_[i] = i;
        fisher_yates(pi_, n);
        for (int i = 0; i < n; i++) invpi_[pi_[i]] = i;
        long long val = compute_I(pi_, n, invpi_);
        if (val > worst) { worst = val; worst_kind = 2; }
    }

    for (int t = 0; t < 500; t++) {
        int h = rand() % n;
        for (int i = 0; i < n; i++) pi_[i] = ((h - i) % n + n) % n;
        int k = 1 + rand() % 3;
        for (int s = 0; s < k; s++) {
            int x = rand() % n, y = rand() % n;
            int tmp = pi_[x]; pi_[x] = pi_[y]; pi_[y] = tmp;
        }
        for (int i = 0; i < n; i++) invpi_[pi_[i]] = i;
        long long val = compute_I(pi_, n, invpi_);
        if (val > worst) { worst = val; worst_kind = 3; }
    }

    printf("n=%d bound=%lld worst=%lld kind=%d %s\n",
           n, bound, worst, worst_kind, worst <= bound ? "OK" : "FAIL");
    free(pi_); free(invpi_);
}

int main(int argc, char **argv) {
    srand(12345);
    if (argc > 1) {
        for (int i = 1; i < argc; i++) scan_n(atoi(argv[i]), 3000);
        return 0;
    }
    int ns[] = {13, 14, 15, 16, 18, 20, 24, 30, 40, 50, 64, 80, 100, 101, 127, 150, 200};
    for (unsigned i = 0; i < sizeof(ns) / sizeof(ns[0]); i++) scan_n(ns[i], 3000);
    return 0;
}
