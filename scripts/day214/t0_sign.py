# exploratory: sign pattern of c_{lam mu} at t=0, s=1/2 (stretch goal: positivity => full support?)
from ds_pointeval import *
from fractions import Fraction as Fr
import random, sys
rng = random.Random(1)
for n in range(2, int(sys.argv[1])+1):
    for lam in parts(n):
        out = coeffs(lam, Fr(1,2), Fr(0), rng)
        full = set(out) == {mu for mu in parts(n) if dominates(mu, lam)}
        print(lam, 'allpos', all(v > 0 for v in out.values()), 'full', full, {mu: str(v) for mu, v in out.items()} if n <= 4 else '', flush=True)
