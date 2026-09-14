# check_C37 — независимая проверка C37–C40

Версия check_C37-1.0, 54.5 с, результат: PASS.

Команда: `python3 checks/check_C37.py --na 6 --nc 7`.

```json
{
 "A": {
  "permutations_and_samples_checked": 944
 },
 "B": {
  "K_n_exact_min_edits": {
   "4": 2,
   "5": 4,
   "6": 6,
   "7": 9,
   "8": 12,
   "9": 16,
   "10": 20
  },
  "all_graphs_max_greedy_edits": {
   "3": 1,
   "4": 2,
   "5": 4,
   "6": 6
  },
  "all_graphs_target": {
   "3": 1,
   "4": 2,
   "5": 4,
   "6": 6
  }
 },
 "C": {
  "equality_cases_are_the_n_reflections": {
   "4": 4,
   "5": 5,
   "6": 6,
   "7": 7,
   "8": 8
  },
  "theorem_3_3": {
   "4": {
    "global_bound": 3,
    "max_I": 2,
    "max_per_pi_bound": 3,
    "target": 2
   },
   "5": {
    "global_bound": 6,
    "max_I": 4,
    "max_per_pi_bound": 6,
    "target": 4
   },
   "6": {
    "global_bound": 9,
    "max_I": 6,
    "max_per_pi_bound": 9,
    "target": 6
   },
   "7": {
    "global_bound": 13,
    "max_I": 9,
    "max_per_pi_bound": 13,
    "target": 9
   },
   "8": {
    "global_bound": 17,
    "max_I": 12,
    "max_per_pi_bound": 17,
    "target": 12
   }
  }
 },
 "D": {
  "reflection_two_cut_value": {
   "4": 2,
   "10": 20,
   "25": 144,
   "40": 380
  },
  "reflections_checked": 37
 },
 "E": {
  "depth1_counterexample": {
   "inv": 21,
   "line": [
    0,
    1,
    8,
    7,
    6,
    5,
    4,
    3,
    2,
    9
   ],
   "n": 10,
   "target": 20,
   "violated_depth1_conditions": 0
  },
  "forall_alpha_exists_beta": {
   "7": {
    "example": [
     0,
     1,
     6,
     5,
     4,
     3,
     2
    ],
    "max_alpha_min_beta": 10,
    "target": 9
   },
   "8": {
    "example": [
     0,
     1,
     7,
     6,
     4,
     5,
     3,
     2
    ],
    "max_alpha_min_beta": 13,
    "target": 12
   }
  },
  "hk_families": {
   "5": {
    "best_family_max": 5,
    "target": 4
   },
   "6": {
    "best_family_max": 7,
    "target": 6
   },
   "7": {
    "best_family_max": 11,
    "target": 9
   },
   "8": {
    "best_family_max": 16,
    "target": 12
   }
  }
 },
 "na": 6,
 "nc": 7,
 "quick": false,
 "result": "PASS",
 "seconds": 54.5,
 "version": "check_C37-1.0"
}
```
