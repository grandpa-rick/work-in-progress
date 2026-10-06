from nonsym import *
from itertools import product
def comps(d, a):
    if a == 1: yield (d,); return
    for i in range(d+1):
        for r in comps(d-i, a-1): yield (i,)+r
ok = True
for a in (2, 3, 4):
    for d in range(0, 4):
        n = a + d
        for beta in comps(d, a):
            gam = tuple(1 + beta[i] - (n if i == 0 else 0) for i in range(a))
            v = ctmono(gam)
            phi = sp.prod([1-t**i for i in range(1, a)])
            pred = (-1)**(a-1)*phi*t**sum(i*beta[i] for i in range(a))
            if sp.simplify(v - pred) != 0: ok = False; print('FAIL', a, beta, sp.factor(v), sp.factor(pred))
print('nonsym 1pt lemma:', ok)
