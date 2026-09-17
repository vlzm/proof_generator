# check_H13I_session9 — новые леммы и расширенный диапазон (H13-I не доказана)

Команда: `python3 checks/check_H13I_session9.py --amax 9`. Версия: check_H13I_session9-1.0.

```text
== check_H13I_session9-1.0 args={'amax': 9, 'amax_inv': 7, 'amax_g': 9}
part 1 n=4: histograms match (max I slow=2 fast=2)
part 1 n=5: histograms match (max I slow=4 fast=4)
part 1 n=6: histograms match (max I slow=6 fast=6)
part 1 n=7: histograms match (max I slow=9 fast=9)
part 1 n=8: histograms match (max I slow=12 fast=12)
part 1 n=9: histograms match (max I slow=16 fast=16)
part 2 (Lemma D, pair-weight identity): PASS (5 random pi for each 3 <= n <= 8, brute sum over all n^2 cuts vs closed form)
part 3 (Lemma E, I(pi)=I(pi^-1)) n=3: 6 pi, 0 mismatches -> PASS
part 3 (Lemma E, I(pi)=I(pi^-1)) n=4: 24 pi, 0 mismatches -> PASS
part 3 (Lemma E, I(pi)=I(pi^-1)) n=5: 120 pi, 0 mismatches -> PASS
part 3 (Lemma E, I(pi)=I(pi^-1)) n=6: 720 pi, 0 mismatches -> PASS
part 3 (Lemma E, I(pi)=I(pi^-1)) n=7: 5040 pi, 0 mismatches -> PASS
part 4 n=4: max G(pi) (position-rotation only) = 2, max J(pi) (value-shift only) = 2, bound = 2 (within / within) -> as recorded
part 4 n=5: max G(pi) (position-rotation only) = 4, max J(pi) (value-shift only) = 4, bound = 4 (within / within) -> as recorded
part 4 n=6: max G(pi) (position-rotation only) = 6, max J(pi) (value-shift only) = 6, bound = 6 (within / within) -> as recorded
part 4 n=7: max G(pi) (position-rotation only) = 10, max J(pi) (value-shift only) = 10, bound = 9 (exceeds / exceeds) -> as recorded
part 4 n=8: max G(pi) (position-rotation only) = 13, max J(pi) (value-shift only) = 13, bound = 12 (exceeds / exceeds) -> as recorded
part 4 n=9: max G(pi) (position-rotation only) = 17, max J(pi) (value-shift only) = 17, bound = 16 (exceeds / exceeds) -> as recorded
part 5 n=4: max I(pi) = 2 (floor((n-1)^2/4) = 2), 24 pi -> ok
part 5 n=5: max I(pi) = 4 (floor((n-1)^2/4) = 4), 120 pi -> ok
part 5 n=6: max I(pi) = 6 (floor((n-1)^2/4) = 6), 720 pi -> ok
part 5 n=7: max I(pi) = 9 (floor((n-1)^2/4) = 9), 5040 pi -> ok
part 5 n=8: max I(pi) = 12 (floor((n-1)^2/4) = 12), 40320 pi -> ok
part 5 n=9: max I(pi) = 16 (floor((n-1)^2/4) = 16), 362880 pi -> ok
part 5 n=10: max I(pi) = 20 (floor((n-1)^2/4) = 20), 3628800 pi -> ok
part 5 n=11: max I(pi) = 25 (floor((n-1)^2/4) = 25), 39916800 pi -> ok
part 5 n=12: max I(pi) = 30 (floor((n-1)^2/4) = 30), 479001600 pi -> ok
part 5 n=13: max I(pi) = 36 (floor((n-1)^2/4) = 36), 6227020800 pi -> ok
total time: 25.9 s
verdict: PASS  (note: this PASS certifies the session-9 lemmas and computational data, NOT H13-I itself, which remains unproven -- see h13_line_model.md par 7)
```
