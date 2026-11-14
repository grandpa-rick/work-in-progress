#!/usr/bin/env python3
"""Recheck the first k=5 positive basis via positivity.negatives (struct() path), and print its full table."""
from positivity import setup, negatives, struct
d, g, nz, r = setup(5)
B = {0: ((),), 1: ((1,),), 2: ((2,), (1, 1)), 3: ((3,), (1, 1, 1)), 4: ((4,), (3, 1), (1, 1, 1, 1)),
     5: ((5,), (2, 2, 1), (1, 1, 1, 1, 1)), 6: ((6,), (5, 1), (3, 3)), 7: ((6, 1), (4, 3)), 8: ((6, 2), (4, 4)),
     9: ((6, 3),), 10: ((6, 4),)}
print('negatives:', negatives(g, B, d))
den = set()
for n in range(2, d + 1):
    for i in range(1, n // 2 + 1):
        for (l, m), c in struct(g, B, i, n - i).items():
            den |= {x.denominator for x in c}
            print(l, '*', m, '=', ' + '.join(f'{x}*s{nu}' for x, nu in zip(c, B[n]) if x))
print('denominators occurring:', sorted(den))
