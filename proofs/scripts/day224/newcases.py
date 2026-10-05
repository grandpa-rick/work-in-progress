"""Genuinely-new l=3, kappa=1 pairs: mu padded to 3 parts has 3 distinct parts and mu_3 = 0 (canonical: no column strip),
and lambda_3 >= 1. Lead = [(s-1)^2] coefficient (v=2 by Thm C). Interpolated from numeric t. Grade: computed."""
from lead2 import *
from kappa import kappa
import sys, time, pickle
nmax = int(sys.argv[1]); out = {}
for n in range(4, nmax+1):
    for lam in parts(n):
        if len(lam) != 3: continue
        mus = [mu for mu in parts(n) if len(mu) == 2 and mu[0] != mu[1] and kappa(lam, mu) == 1]
        if not mus: continue
        cache = {}
        def val(mu, t):
            if t not in cache: cache[t] = second(*lam, t)
            return cache[t].get(mu, 0)
        for mu in mus:
            p = interp(lambda t: val(mu, t), npts=3*n+12)
            out[(lam, mu)] = p
            print(lam, '->', mu, ':', p, ' | t=0:', p.subs(ts, 0), ' t=1:', p.subs(ts, 1), flush=True)
pickle.dump(out, open(f'newcases_n{nmax}.pkl', 'wb'))
