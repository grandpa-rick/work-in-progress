from lead2 import *
import sys
nmax = int(sys.argv[1])
for n in range(3, nmax+1):
    for a in range(1, n-1):
        for b in range(1, n-a):
            c = n-a-b
            if c > b: continue
            cache = {}
            def g(t):
                if t not in cache: cache[t] = toe(Gam(a, b, c, t), t)
                return cache[t]
            keys = set()
            for tt in (2, 3, 5): keys |= set(g(Fr(tt)).keys())
            res = {mu: interp(lambda t: g(t).get(mu, 0), 2*n*n+8) for mu in sorted(keys, reverse=True)}
            print(f'Gam_{a}(e{b},e{c}):', {k: v for k, v in res.items() if v != 0}, flush=True)
