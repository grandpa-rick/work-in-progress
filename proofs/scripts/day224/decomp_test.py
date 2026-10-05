"""Conjecture (v=2 multiplicativity): for kappa(lam,mu)=l-2, Lead = sum over decompositions realizing kappa of prod of block leads."""
from second_gen import *
from kappa import kappa, dom
from collections import Counter
import itertools, sys
t = Fr(int(sys.argv[2])) if len(sys.argv) > 2 else Fr(3)
def setparts(lst):
    if not lst: yield []; return
    f = lst[0]
    for p in setparts(lst[1:]):
        for i in range(len(p)): yield p[:i] + [[f]+p[i]] + p[i+1:]
        yield [[f]] + p
def msplits(mu, sizes):
    """ordered splits of multiset mu into sub-multisets with given sums."""
    if not sizes:
        if not mu: yield []
        return
    mu = list(mu)
    seen = set()
    for r in range(1, len(mu)+1):
        for comb in itertools.combinations(range(len(mu)), r):
            sub = tuple(sorted([mu[i] for i in comb], reverse=True))
            if sum(sub) != sizes[0] or sub in seen: continue
            seen.add(sub)
            rest = [mu[i] for i in range(len(mu)) if i not in comb]
            for tail in msplits(rest, sizes[1:]): yield [sub] + tail
_lead_cache = {}
def block_lead(lamC, nu):
    lamC = tuple(sorted(lamC, reverse=True))
    if len(lamC) == 1: return Fr(1) if nu == lamC else None
    if not dom(nu, lamC) or nu == lamC: return None
    if len(lamC) == 2:
        if kappa(lamC, nu) != 1: return None
        return Mthm(lamC[0], lamC[1], t).get(nu, Fr(0))   # (s-1)^1 coefficient = lead at v=1
    if len(lamC) == 3:
        if kappa(lamC, nu) != 1: return None
        k = (lamC, nu)
        if k not in _lead_cache: _lead_cache[k] = second(*lamC, t).get(nu, Fr(0))
        return _lead_cache[k]
    return None
nmax = int(sys.argv[1]); ok = bad = 0
for n in range(4, nmax+1):
    for lam in parts(n):
        l = len(lam)
        if l < 4: continue
        S2 = second_general(lam, t)[2]
        for mu in parts(n):
            if kappa(lam, mu) != l-2: continue
            pred = Fr(0)
            for pi in setparts(list(range(l))):
                if len(pi) != l-2: continue
                sizes = [sum(lam[i] for i in C) for C in pi]
                for sp in msplits(mu, sizes):
                    prod = Fr(1)
                    for C, nu in zip(pi, sp):
                        bl = block_lead([lam[i] for i in C], nu)
                        if bl is None: prod = None; break
                        prod *= bl
                    if prod is not None: pred += prod
            got = S2.get(mu, Fr(0))
            g = (got == pred); ok += g; bad += not g
            if not g: print('FAIL', lam, mu, got, pred, flush=True)
print('t=', t, 'ok', ok, 'bad', bad)
