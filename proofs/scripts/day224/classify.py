"""Classify l=3, kappa=1, two-part mu (distinct parts) pairs by which PROVED formula covers them."""
from hl import parts
from kappa import kappa
cnt = {}
for n in range(4, 13):
    for lam in parts(n):
        if len(lam) != 3: continue
        for mu in parts(n):
            if len(mu) != 2 or mu[0] == mu[1] or kappa(lam, mu) != 1: continue
            x, y = mu
            if y == 1 or x - y == 1: cls = 'N1 (mu~(n-1,1))'
            elif lam[2] == 1 or lam[0] == x - 1: cls = 'lambda-has-1 (a=1 formula)'
            else: cls = 'OPEN'
            cnt[cls] = cnt.get(cls, 0) + 1
            if cls == 'OPEN' and n <= 11: print(n, lam, mu)
print(cnt)
