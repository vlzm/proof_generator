# check_C28 — PASS

Версии: check_C28-1.0-pr5, strict_upper-1.0, oracle-1.0; команда `checks/check_C28_pr5.py`; 19.0 с.

Покрытие: reflections_all_h: 4<=n<=64, all h (2074 inputs), chosen shift, words built, executed and reduced; min_over_all_c: 4<=n<=40, all h, all c (report only); sigma_n: 4<=n<=200

Все assert'ы лемм 2–6, 8 и случаев 1–6 прошли. Значения `len − B_n` по случаям (выбранный сдвиг):

- случай 1: [-1, 0]
- случай 2: [-1]
- случай 3: [0]
- случай 4: [0]
- случай 5: [1]
- случай 6: [-1]

Максимум `len − B_n` = 1, максимум `red − B_n` = 0.
Минимум по всем c при n <= 40: max `min_c len_R − B_n` = 1, max `min_c red_R − B_n` = 0; выбранный сдвиг не минимален по длине на 395 входах (это не влияет на теорему).
sigma_n при 4 <= n <= 200: `len − B_n` принимает значения [0], N_X = floor((n−1)^2/4) на всех n.
