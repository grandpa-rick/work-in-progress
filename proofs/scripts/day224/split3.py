from lead2 import *
from kappa import kappa
import sys, itertools
nmax = int(sys.argv[1])
for n in range(4, nmax+1):
    for lam in parts(n):
        if len(lam) != 3: continue
        mus = [mu for mu in parts(n) if mu != (n,) and kappa(lam, mu) == 1]
        if not mus: continue
        for perm in sorted(set(itertools.permutations(lam))):
            cache = {}
            def pc(t):
                if t not in cache: cache[t] = pieces(*perm, t)
                return cache[t]
            for mu in mus:
                G = interp(lambda t: pc(t)['Gam'].get(mu, 0), 3*n+12)
                DD = interp(lambda t: pc(t)['DD'].get(mu, 0), 3*n+12)
                other = [k for k in ('Ebc', 'Ecb', 'eE') if interp(lambda t: pc(t)[k].get(mu, 0), 3*n+12) != 0]
                print(perm, '->', mu, ' Gam:', G, ' | DD:', DD, ' | others nonzero:', other, flush=True)
