# H13-I session 9: quadrant reformulation and rejected attempts (h13_i_quadrant-1.0)

nmax = 8, elapsed = 37.1 s

## Part 1: sum_{a,b} concordant(a,b) = sum_{pairs} K(g,h) -- identity, PROVED

- n=4: 24 perms, identity holds: True
- n=5: 120 perms, identity holds: True
- n=6: 720 perms, identity holds: True
- n=7: 5040 perms, identity holds: True
- n=8: 40320 perms, identity holds: True

## Part 2: second-moment bound S2/S -- REFUTED (fails already at the equality case)

- n=4: target=4, worst pi=[0, 3, 2, 1], S2/S=3.400, actual max=4, bound holds: False
- n=5: target=6, worst pi=[0, 4, 3, 2, 1], S2/S=5.200, actual max=6, bound holds: False
- n=6: target=9, worst pi=[0, 5, 4, 3, 2, 1], S2/S=7.400, actual max=9, bound holds: False
- n=7: target=12, worst pi=[0, 6, 5, 4, 3, 2, 1], S2/S=10.000, actual max=12, bound holds: False
- n=8: target=16, worst pi=[0, 7, 6, 5, 4, 3, 2, 1], S2/S=13.000, actual max=16, bound holds: False

## Part 3: joint (row+column) fixed points -- REFUTED as sufficient

- n=4: 0/24 perms have a joint fixed point exceeding floor((n-1)^2/4)=2
- n=5: 0/120 perms have a joint fixed point exceeding floor((n-1)^2/4)=4
- n=6: 0/720 perms have a joint fixed point exceeding floor((n-1)^2/4)=6
- n=7: 49/5040 perms have a joint fixed point exceeding floor((n-1)^2/4)=9
  example: pi=[0, 1, 6, 5, 4, 3, 2], worst fixed point inv=10, all fixed points=[(0, 5, 8), (0, 6, 8), (1, 1, 10), (2, 3, 8), (2, 4, 8), (3, 2, 8), (3, 3, 8), (4, 2, 8), (5, 0, 8), (6, 0, 8), (6, 6, 8)]
- n=8: 0/40320 perms have a joint fixed point exceeding floor((n-1)^2/4)=12

## Part 4: 321-avoidance of optimal cuts -- REFUTED

- n=4: 0/24 perms have no 321-avoiding optimal cut
- n=5: 5/120 perms have no 321-avoiding optimal cut
  example: pi=[0, 4, 3, 2, 1], best inv=4, an optimal line=[2, 1, 0, 4, 3]
- n=6: 150/720 perms have no 321-avoiding optimal cut
  example: pi=[0, 1, 2, 5, 4, 3], best inv=3, an optimal line=[0, 1, 2, 5, 4, 3]
- n=7: 1785/5040 perms have no 321-avoiding optimal cut
  example: pi=[0, 1, 2, 3, 6, 5, 4], best inv=3, an optimal line=[0, 1, 2, 3, 6, 5, 4]
- n=8: 24152/40320 perms have no 321-avoiding optimal cut
  example: pi=[0, 1, 2, 3, 4, 7, 6, 5], best inv=3, an optimal line=[0, 1, 2, 3, 4, 7, 6, 5]
