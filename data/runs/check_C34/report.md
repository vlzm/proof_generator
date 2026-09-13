# check_C34 — линейная модель: универсальный маршрут, антиподальная транспозиция, профиль

Команда: `python3 checks/check_C34.py --amax 8`. Версии: check_C34-1.0, line_cocktail-1.0, oracle-1.0.

```text
== check_C34-1.0 line_cocktail-1.0 oracle-1.0 args={'amax': 8}
part A n=4: 24 lines sorted, N_X = inv, N_rot = B_n - 1 = 5, 0.0 s
part A n=5: 120 lines sorted, N_X = inv, N_rot = B_n - 1 = 9, 0.0 s
part A n=6: 720 lines sorted, N_X = inv, N_rot = B_n - 1 = 14, 0.0 s
part A n=7: 5040 lines sorted, N_X = inv, N_rot = B_n - 1 = 20, 0.1 s
part A n=8: 40320 lines sorted, N_X = inv, N_rot = B_n - 1 = 27, 1.3 s
part B n=4: tau=(0, 1, 3, 2) d=5 (n+1=5) I=1 word LLXRR sorts=True -> ok
part B n=6: tau=(0, 1, 2, 4, 3, 5) d=7 (n+1=7) I=1 word LLLXRRR sorts=True -> ok
part B n=8: tau=(0, 1, 2, 3, 5, 4, 6, 7) d=9 (n+1=9) I=1 word LLLLXRRRR sorts=True -> ok
part B n=10: tau=(0, 1, 2, 3, 4, 6, 5, 7, 8, 9) d=11 (n+1=11) I=1 word LLLLLXRRRRR sorts=True -> ok
part C n=4: max(d - I) = 4 (floor(n^2/4) = 4), max I = 2 (floor((n-1)^2/4) = 2), max(d - 2I) = 3 -> ok
part C n=5: max(d - I) = 6 (floor(n^2/4) = 6), max I = 4 (floor((n-1)^2/4) = 4), max(d - 2I) = 3 -> ok
part C n=6: max(d - I) = 9 (floor(n^2/4) = 9), max I = 6 (floor((n-1)^2/4) = 6), max(d - 2I) = 5 -> ok
part C n=7: max(d - I) = 12 (floor(n^2/4) = 12), max I = 9 (floor((n-1)^2/4) = 9), max(d - 2I) = 6 -> ok
part C n=8: max(d - I) = 16 (floor(n^2/4) = 16), max I = 12 (floor((n-1)^2/4) = 12), max(d - 2I) = 8 -> ok
part C n=9: max(d - I) = 20 (floor(n^2/4) = 20), max I = 16 (floor((n-1)^2/4) = 16), max(d - 2I) = 10 -> ok
part C n=10: max(d - I) = 25 (floor(n^2/4) = 25), max I = 20 (floor((n-1)^2/4) = 20), max(d - 2I) = 13 -> ok
verdict: PASS
```
