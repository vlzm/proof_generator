# check_H13I — PASS

Версия check_H13I-1.0; brute-force диапазон 4 <= n <= 8; время 49.4 с.

- A: shift recurrences verified on 5713488 cut pairs, 4 <= n <= 8: PASS
- B: duality I(pi) + M(pi^r) = C(n,2) on 46224 permutations (4 <= n <= 8); arithmetic identity 4 <= n <= 200: PASS
- C: winding lemma on 46224 permutations, 4 <= n <= 8: PASS
- D: reflections I = floor((n-1)^2/4) (brute force 4 <= n <= 9, closed form 4 <= n <= 200): PASS
- E: run decomposition at a descent on 14652 (pi, cut) pairs, 4 <= n <= 6: PASS
- F1: brute force H13-I with equality exactly on reflections, 4 <= n <= 8: PASS
- F2: toric_inv-1.0 JSONs for n in [4, 5, 6, 7, 8, 9, 10, 11, 12, 13]: PASS
- F3: session-8 line_profile-1.0 histograms agree for n in [4, 5, 6, 7, 8, 9, 10]: PASS
- G: pair-area and average-concordance formulas on 864 permutations, 4 <= n <= 6: PASS
