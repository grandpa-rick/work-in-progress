#!/usr/bin/env python3
"""Exhaustive search for Schur-subset bases (one per degree) of G_(2^k) with ALL structure constants >= 0.
Same math as positivity.search, with caching: nonneg[(B_n, lam, mu)] computed once.
Also reports, per degree, how many B_n survive for SOME choice of lower bases.  Grade: computed."""
import sys, itertools, time
from engine import *
from positivity import setup, bases_in_degree

import os
LIMIT = int(os.environ.get("LIMIT", 10**6))
for k in map(int, sys.argv[1:]):
    t0 = time.time()
    d, g, nz, r = setup(k)
    choices = {n: bases_in_degree(g, nz, r, n) for n in range(d + 1)}
    cache = {}
    def nonneg(Bn, n, lam, mu):
        key = (Bn, lam, mu)
        if key not in cache:
            x = solve([g.img(l) for l in Bn], g.img_of(lr(lam, mu), n))
            cache[key] = all(t >= 0 for t in x)
        return cache[key]
    # backtracking with memo on "state" = (B_1..B_{n-1}) is exponential; prune by degree.
    sols = [0]; first = []
    surv = {n: set() for n in range(d + 1)}
    def rec(n, B):
        if sols[0] >= LIMIT: return
        if n > d:
            sols[0] += 1
            if len(first) < 3: first.append(dict(B))
            return
        for S in choices[n]:
            good = all(nonneg(S, n, lam, mu) for i in range(1, n // 2 + 1) for lam in B[i] for mu in B[n - i])
            if good:
                surv[n].add(S); B[n] = S; rec(n + 1, B); del B[n]
    rec(0, {})
    print('=' * 60)
    print('k =', k, 'Hilb', [r[n] for n in range(d + 1)])
    print(' #Schur bases per degree:', [len(choices[n]) for n in range(d + 1)])
    print(' #bases per degree appearing in some positive partial chain:', [len(surv[n]) for n in range(d + 1)])
    print(" # fully POSITIVE Schur-subset bases:", sols[0], ' time %.1fs' % (time.time() - t0))
    for s in first: print('   e.g.', {n: s[n] for n in range(2, d + 1)})
    sys.stdout.flush()
