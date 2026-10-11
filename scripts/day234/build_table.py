# Build ground-truth table c[lam][mu] for |lam|<=7 via engine.py (subset formula only), N=7 variables.
# Also stability control: recompute n<=6 with N=n and N=n+2 and compare.
import sys, pickle, time; sys.path.insert(0, '.')
from engine import *
T = {}
E7 = Engine(7, seed=7)
t0 = time.time()
for n in range(1, 8):
    for lam in partitions(n):
        T[lam] = E7.estar(lam)
print('table n<=7 built', round(time.time()-t0,1), 's', sum(1 for _ in T))
pickle.dump({lam: {mu: str(c.as_expr()) for mu, c in d.items()} for lam, d in T.items()}, open('ctable_n7.pkl', 'wb'))
bad = 0; cnt = 0
for n in range(1, 7):
    for N in (n, n + 2):
        En = Engine(N, seed=100 + N)
        for lam in partitions(n):
            G = En.estar(lam); cnt += 1
            if G != T[lam]: bad += 1; print('STABILITY MISMATCH', N, lam)
print('N-stability (N=n, n+2 vs 7), n<=6:', cnt, 'comparisons,', bad, 'mismatches')
