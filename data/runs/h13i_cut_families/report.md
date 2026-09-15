# h13i_cut_families — ограниченные семейства разрезов и индукция

Цель: проверить, достаточно ли для H13-I какого-нибудь якорного семейства из n разрезов, и проходит ли индукция удалением k элементов. Покрытие: исчерпывающе по всем pi при 4 <= n <= 8 (часть 1, 3) и 5 <= n <= 7 (часть 2), плюс именованные аффинные свидетели при n = 8, 9. Seed не используется (перебор).

Команда: `python3 experiments/h13i_cut_families.py --nmax 7 --nmax-anchor 8`. Версии: h13i_cut_families-1.0, oracle-1.0.

```text
== h13i_cut_families-1.0 oracle-1.0 args={'nmax': 7, 'nmax_anchor': 8}
part 1 n=4 target=2: anchored families reaching the target: 0 of 16; best family (0, 0) gives 3 on pi=(0, 3, 2, 1), 0.0 s
part 1 n=5 target=4: anchored families reaching the target: 0 of 25; best family (0, 1) gives 5 on pi=(0, 3, 1, 4, 2), 0.0 s
part 1 n=6 target=6: anchored families reaching the target: 0 of 36; best family (0, 1) gives 7 on pi=(0, 5, 4, 3, 2, 1), 0.0 s
part 1 n=7 target=9: anchored families reaching the target: 0 of 49; best family (0, 1) gives 11 on pi=(0, 5, 3, 1, 6, 4, 2), 0.3 s
part 1 n=8 target=12: anchored families reaching the target: 0 of 64; best family (0, 1) gives 16 on pi=(0, 7, 6, 5, 4, 3, 2, 1), 3.5 s
part 2 n=5 delete 1: allowed gap 2, worst needed 2 on pi=(0, 2, 4, 1, 3) -> ok, 0.0 s
part 2 n=6 delete 1: allowed gap 2, worst needed 2 on pi=(0, 1, 4, 5, 2, 3) -> ok, 0.0 s
part 2 n=7 delete 1: allowed gap 3, worst needed 3 on pi=(0, 2, 4, 6, 1, 3, 5) -> ok, 0.5 s
part 2 n=6 delete 2: allowed gap 4, worst needed 4 on pi=(0, 5, 4, 3, 2, 1) -> ok, 0.1 s
part 2 n=7 delete 2: allowed gap 5, worst needed 5 on pi=(0, 6, 5, 4, 3, 2, 1) -> ok, 1.0 s
part 2 n=7 delete 3: allowed gap 7, worst needed 7 on pi=(0, 6, 5, 4, 3, 2, 1) -> ok, 1.4 s
part 2 witness (0, 5, 2, 7, 4, 1, 6, 3): I=11, delete 1: best I(pi')=6 at (7,), allowed gap 3, needed 5 -> FAILS
part 2 witness (0, 5, 2, 7, 4, 1, 6, 3): I=11, delete 2: best I(pi')=4 at (6, 7), allowed gap 6, needed 7 -> FAILS
part 2 witness (0, 5, 2, 7, 4, 1, 6, 3): I=11, delete 3: best I(pi')=3 at (5, 6, 7), allowed gap 8, needed 8 -> ok
part 2 witness (0, 3, 6, 1, 4, 7, 2, 5): I=9, delete 1: best I(pi')=5 at (7,), allowed gap 3, needed 4 -> FAILS
part 2 witness (0, 3, 6, 1, 4, 7, 2, 5): I=9, delete 2: best I(pi')=4 at (5, 7), allowed gap 6, needed 5 -> ok
part 2 witness (0, 3, 6, 1, 4, 7, 2, 5): I=9, delete 3: best I(pi')=3 at (5, 6, 7), allowed gap 8, needed 6 -> ok
part 2 witness (0, 2, 4, 6, 8, 1, 3, 5, 7): I=10, delete 1: best I(pi')=6 at (8,), allowed gap 4, needed 4 -> ok
part 2 witness (0, 2, 4, 6, 8, 1, 3, 5, 7): I=10, delete 2: best I(pi')=6 at (4, 8), allowed gap 7, needed 4 -> ok
part 2 witness (0, 2, 4, 6, 8, 1, 3, 5, 7): I=10, delete 3: best I(pi')=4 at (2, 5, 8), allowed gap 10, needed 6 -> ok
part 2 witness (0, 5, 1, 6, 2, 7, 3, 8, 4): I=10, delete 1: best I(pi')=6 at (8,), allowed gap 4, needed 4 -> ok
part 2 witness (0, 5, 1, 6, 2, 7, 3, 8, 4): I=10, delete 2: best I(pi')=6 at (7, 8), allowed gap 7, needed 4 -> ok
part 2 witness (0, 5, 1, 6, 2, 7, 3, 8, 4): I=10, delete 3: best I(pi')=4 at (2, 5, 8), allowed gap 10, needed 6 -> ok
part 3 n=4: min over pi of #cuts with inv <= 2 is 4 (of 16) on pi=(0, 1, 2, 3), 0.0 s
part 3 n=5: min over pi of #cuts with inv <= 4 is 7 (of 25) on pi=(0, 1, 4, 3, 2), 0.0 s
part 3 n=6: min over pi of #cuts with inv <= 6 is 6 (of 36) on pi=(0, 5, 4, 3, 2, 1), 0.0 s
part 3 n=7: min over pi of #cuts with inv <= 9 is 10 (of 49) on pi=(0, 1, 6, 5, 4, 3, 2), 0.1 s
part 3 n=8: min over pi of #cuts with inv <= 12 is 8 (of 64) on pi=(0, 7, 6, 5, 4, 3, 2, 1), 0.6 s
total 7.5 s, peak 12 MB
```
