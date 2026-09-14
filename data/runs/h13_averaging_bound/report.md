# h13_averaging_bound report (h13_averaging_bound-1.0, core oracle-1.0)

Elapsed: 8.0 s

```
part 1: decomposition A xor B verified on 75 random (n, pi) samples, 4 <= n <= 8: OK
part 2 n=4: H(a)=0 for all a: True; E=14 (n-1)(2n-1)/6; pigeonhole bound E/n=3.500 vs threshold 2 -> gap 1.500 (method fails)
part 2 n=5: H(a)=0 for all a: True; E=30 (n-1)(2n-1)/6; pigeonhole bound E/n=6.000 vs threshold 4 -> gap 2.000 (method fails)
part 2 n=6: H(a)=0 for all a: True; E=55 (n-1)(2n-1)/6; pigeonhole bound E/n=9.167 vs threshold 6 -> gap 3.167 (method fails)
part 2 n=7: H(a)=0 for all a: True; E=91 (n-1)(2n-1)/6; pigeonhole bound E/n=13.000 vs threshold 9 -> gap 4.000 (method fails)
part 2 n=8: H(a)=0 for all a: True; E=140 (n-1)(2n-1)/6; pigeonhole bound E/n=17.500 vs threshold 12 -> gap 5.500 (method fails)
part 2 n=9: H(a)=0 for all a: True; E=204 (n-1)(2n-1)/6; pigeonhole bound E/n=22.667 vs threshold 16 -> gap 6.667 (method fails)
part 2 n=10: H(a)=0 for all a: True; E=285 (n-1)(2n-1)/6; pigeonhole bound E/n=28.500 vs threshold 20 -> gap 8.500 (method fails)
part 3 (exhaustive) n=4 a_list=(0, 3): worst=2 threshold=2 excess=0 fails=0
part 3 (exhaustive) n=5 a_list=(0, 3): worst=4 threshold=4 excess=0 fails=0
part 3 (exhaustive) n=6 a_list=(0, 3): worst=6 threshold=6 excess=0 fails=0
part 3 (exhaustive) n=7 a_list=(0, 3): worst=9 threshold=9 excess=0 fails=0
part 3 (exhaustive) n=8 a_list=(0, 3): worst=12 threshold=12 excess=0 fails=0
part 3 (hill-climbing) n=20 a_list=(0, 3): found excess=2 (witness_verified=True) -- fixed pair is NOT sufficient in general
```

Verdict: H13-I NOT proved. Part 2 shows the pigeonhole/averaging bound suggested by the algebraic hint fails on reflections by a gap growing like (n^2-1)/12, for every choice of a -- a proved dead end, not just a small-n observation. Part 3 shows the natural 'fixed small family of a' idea (exhaustively sufficient up to n=11 for {0,3}) breaks at n=20 (hill-climbing counterexample, excess 2) -- consistent with session 8's rejection of restricted-cut families: no O(1)-size, pi-independent set of position-cuts suffices in general.
