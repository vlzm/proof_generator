# h13_attempts report (session 9)

```
== h13_attempts-1.0 oracle-1.0 args={'amax': 8, 'rmax': 12}
route1 n=4: bound=2, max over pi of value-rotation-only min = 2 (pi=(0, 2, 3, 1)), max over pi of position-rotation-only min = 2 (pi=(0, 3, 1, 2)) -> both within bound
route1 n=5: bound=4, max over pi of value-rotation-only min = 4 (pi=(0, 2, 4, 3, 1)), max over pi of position-rotation-only min = 4 (pi=(0, 4, 1, 3, 2)) -> both within bound
route1 n=6: bound=6, max over pi of value-rotation-only min = 6 (pi=(0, 2, 4, 5, 3, 1)), max over pi of position-rotation-only min = 6 (pi=(0, 4, 3, 2, 1, 5)) -> both within bound
route1 n=7: bound=9, max over pi of value-rotation-only min = 10 (pi=(0, 5, 4, 3, 2, 1, 6)), max over pi of position-rotation-only min = 10 (pi=(0, 5, 4, 3, 2, 1, 6)) -> EXCEEDS bound (route fails)
route1 n=8: bound=12, max over pi of value-rotation-only min = 13 (pi=(0, 6, 5, 3, 4, 2, 1, 7)), max over pi of position-rotation-only min = 13 (pi=(0, 6, 5, 3, 4, 2, 1, 7)) -> EXCEEDS bound (route fails)
route2 n=4: bound=2, true I(rev_n)=2, mean=3.50, max=6, Cantelli guaranteed upper bound=2.60 -> NOT tight enough
route2 n=5: bound=4, true I(rev_n)=4, mean=6.00, max=10, Cantelli guaranteed upper bound=4.80 -> NOT tight enough
route2 n=6: bound=6, true I(rev_n)=6, mean=9.17, max=15, Cantelli guaranteed upper bound=7.60 -> NOT tight enough
route2 n=7: bound=9, true I(rev_n)=9, mean=13.00, max=21, Cantelli guaranteed upper bound=11.00 -> NOT tight enough
route2 n=8: bound=12, true I(rev_n)=12, mean=17.50, max=28, Cantelli guaranteed upper bound=15.00 -> NOT tight enough
route2 n=9: bound=16, true I(rev_n)=16, mean=22.67, max=36, Cantelli guaranteed upper bound=19.60 -> NOT tight enough
route2 n=10: bound=20, true I(rev_n)=20, mean=28.50, max=45, Cantelli guaranteed upper bound=24.80 -> NOT tight enough
route2 n=11: bound=25, true I(rev_n)=25, mean=35.00, max=55, Cantelli guaranteed upper bound=30.60 -> NOT tight enough
route2 n=12: bound=30, true I(rev_n)=30, mean=42.17, max=66, Cantelli guaranteed upper bound=37.00 -> NOT tight enough
route1 sufficed for all n checked: False (expected False -- both single-axis restrictions fail already at n=7)
total time: 1.9 s
```
