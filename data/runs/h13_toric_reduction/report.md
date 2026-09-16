# h13_toric_reduction — нестед-волк тождество и провал редуцированных семейств разрезов

Команда: `python3 experiments/h13_toric_reduction.py --nmax 7 --identity-nmax 20`. Версия: h13_toric_reduction-1.0, oracle-1.0.

```text
== h13_toric_reduction-1.0 oracle-1.0 args={'nmax': 7, 'identity_nmax': 20, 'identity_trials': 80}
identities n=4: (W-A) and (W-B) hold exactly on 80 random perms, all a, all k
identities n=5: (W-A) and (W-B) hold exactly on 80 random perms, all a, all k
identities n=6: (W-A) and (W-B) hold exactly on 80 random perms, all a, all k
identities n=7: (W-A) and (W-B) hold exactly on 80 random perms, all a, all k
identities n=8: (W-A) and (W-B) hold exactly on 80 random perms, all a, all k
identities n=9: (W-A) and (W-B) hold exactly on 80 random perms, all a, all k
identities n=10: (W-A) and (W-B) hold exactly on 80 random perms, all a, all k
identities n=11: (W-A) and (W-B) hold exactly on 80 random perms, all a, all k
identities n=12: (W-A) and (W-B) hold exactly on 80 random perms, all a, all k
identities n=13: (W-A) and (W-B) hold exactly on 80 random perms, all a, all k
identities n=14: (W-A) and (W-B) hold exactly on 80 random perms, all a, all k
identities n=15: (W-A) and (W-B) hold exactly on 80 random perms, all a, all k
identities n=16: (W-A) and (W-B) hold exactly on 80 random perms, all a, all k
identities n=17: (W-A) and (W-B) hold exactly on 80 random perms, all a, all k
identities n=18: (W-A) and (W-B) hold exactly on 80 random perms, all a, all k
identities n=19: (W-A) and (W-B) hold exactly on 80 random perms, all a, all k
identities n=20: (W-A) and (W-B) hold exactly on 80 random perms, all a, all k
family [a=0, any b] n=4: worst(min over family) = 2, floor((n-1)^2/4) = 2 -> OK (0.0 s)
family [a=0, any b] n=5: worst(min over family) = 4, floor((n-1)^2/4) = 4 -> OK (0.0 s)
family [a=0, any b] n=6: worst(min over family) = 6, floor((n-1)^2/4) = 6 -> OK (0.0 s)
family [a=0, any b] n=7: worst(min over family) = 10, floor((n-1)^2/4) = 9 -> FAIL, counterexample (0, 5, 4, 3, 2, 1, 6) (0.1 s)
family [a=0, any b]: first fails at n=7
family [b=0, any a] n=4: worst(min over family) = 2, floor((n-1)^2/4) = 2 -> OK (0.0 s)
family [b=0, any a] n=5: worst(min over family) = 4, floor((n-1)^2/4) = 4 -> OK (0.0 s)
family [b=0, any a] n=6: worst(min over family) = 6, floor((n-1)^2/4) = 6 -> OK (0.0 s)
family [b=0, any a] n=7: worst(min over family) = 10, floor((n-1)^2/4) = 9 -> FAIL, counterexample (0, 5, 4, 3, 2, 1, 6) (0.1 s)
family [b=0, any a]: first fails at n=7
family [b=a] n=4: worst(min over family) = 4, floor((n-1)^2/4) = 2 -> FAIL, counterexample (2, 3, 0, 1) (0.0 s)
family [b=a] n=5: worst(min over family) = 6, floor((n-1)^2/4) = 4 -> FAIL, counterexample (2, 3, 4, 0, 1) (0.0 s)
family [b=a] n=6: worst(min over family) = 10, floor((n-1)^2/4) = 6 -> FAIL, counterexample (2, 5, 4, 1, 0, 3) (0.0 s)
family [b=a] n=7: worst(min over family) = 12, floor((n-1)^2/4) = 9 -> FAIL, counterexample (2, 3, 6, 5, 1, 0, 4) (0.1 s)
family [b=a]: first fails at n=4
family [b=-a] n=4: worst(min over family) = 6, floor((n-1)^2/4) = 2 -> FAIL, counterexample (3, 2, 1, 0) (0.0 s)
family [b=-a] n=5: worst(min over family) = 10, floor((n-1)^2/4) = 4 -> FAIL, counterexample (4, 3, 2, 1, 0) (0.0 s)
family [b=-a] n=6: worst(min over family) = 15, floor((n-1)^2/4) = 6 -> FAIL, counterexample (5, 4, 3, 2, 1, 0) (0.0 s)
family [b=-a] n=7: worst(min over family) = 21, floor((n-1)^2/4) = 9 -> FAIL, counterexample (6, 5, 4, 3, 2, 1, 0) (0.1 s)
family [b=-a]: first fails at n=4
family [union(a=0,any b ; b=0,any a)] n=4: worst(min over family) = 2, floor((n-1)^2/4) = 2 -> OK (0.0 s)
family [union(a=0,any b ; b=0,any a)] n=5: worst(min over family) = 4, floor((n-1)^2/4) = 4 -> OK (0.0 s)
family [union(a=0,any b ; b=0,any a)] n=6: worst(min over family) = 6, floor((n-1)^2/4) = 6 -> OK (0.0 s)
family [union(a=0,any b ; b=0,any a)] n=7: worst(min over family) = 10, floor((n-1)^2/4) = 9 -> FAIL, counterexample (0, 5, 4, 3, 2, 1, 6) (0.2 s)
family [union(a=0,any b ; b=0,any a)]: first fails at n=7
```
