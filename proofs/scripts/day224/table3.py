from lead2 import *
from kappa import kappa
import sys, time
nmax = int(sys.argv[1])
for n in range(3, nmax+1):
    for lam in parts(n):
        if len(lam) != 3: continue
        mus = [mu for mu in parts(n) if mu != (n,) and kappa(lam, mu) == 1]
        if not mus: continue
        t0 = time.time()
        cache = {}
        def val(mu, t):
            if t not in cache: cache[t] = second(*lam, t)
            return cache[t].get(mu, 0)
        for mu in mus:
            p = interp(lambda t: val(mu, t), npts=3*n+12)
            print(lam, '->', mu, ':', p, ' | t=0:', p.subs(ts, 0), ' t=1:', p.subs(ts, 1), flush=True)
