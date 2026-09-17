# check_H13I report (check_H13I-1.0)

part A: 17904 permutations, O(n^2) recurrence vs brute force, PASS
part B1 (value-shift only, a=0 fixed): worst={4: 2, 5: 4, 6: 6, 7: 10, 8: 13}, exceeds bound 9 at n=7 as expected -- rejects this simplification
part B2 (point-aligned cuts b=pi(a)): worst={4: 3, 5: 6, 6: 10, 7: 15, 8: 21}, bounds=[2, 4, 6, 9, 12]
part B3 (diagonal cuts b=a): worst={4: 4, 5: 6, 6: 10, 7: 12, 8: 17}, bounds=[2, 4, 6, 9, 12]
part C: n=11 exhaustive (39916800 permutations, 8s): max_I=25 == bound=25, exceed_count=0 -- PASS
part C: n=12 exhaustive (479001600 permutations, 120s): max_I=30 == bound=30, exceed_count=0 -- PASS
part D: sa_search.log scanned, max margin seen = 0 (<=0 expected), PASS

OVERALL: PASS
