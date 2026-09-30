# exploratory: e_k * e_mu at t=0, s=1/2 : positive on full up-set of mu∪k ?
from ds_pointeval import *
from fractions import Fraction as Fr
import random, math, sys
rng = random.Random(2)
def op(k, mu, x, s, t):
    m = len(x); tot = Fr(0)
    for A in itertools.combinations(range(m), k):
        As = set(A); w = Fr(1)
        for i in A:
            w *= x[i]
            for j in range(m):
                if j not in As: w *= (x[i]-t*x[j])/(x[i]-x[j])
        y = [s*x[i] if i in As else x[i] for i in range(m)]
        tot += w*math.prod([esym(p, y) for p in mu])
    return tot
def opcoeffs(k, mu, s, t):
    n = k+sum(mu); P = list(parts(n))
    pts = [tuple(Fr(rng.randint(-10**6, 10**6), rng.randint(1, 10**4)) for _ in range(n)) for _ in range(len(P))]
    M = [[math.prod([esym(p, x) for p in nu]) for nu in P] for x in pts]
    c = solve(M, [op(k, mu, x, s, t) for x in pts]); return {nu: ci for nu, ci in zip(P, c) if ci != 0}
if __name__ == '__main__':
    bad = 0
    for n in range(2, int(sys.argv[1])+1):
        for k in range(1, n):
            for mu in parts(n-k):
                lam = tuple(sorted(mu+(k,), reverse=True))
                out = opcoeffs(k, mu, Fr(1,2), Fr(0))
                up = {nu for nu in parts(n) if dominates(nu, lam)}
                pos = all(v > 0 for v in out.values()); full = set(out) == up; sub = set(out) <= up
                if not (pos and full): bad += 1; print('k', k, 'mu', mu, 'pos', pos, 'full', full, 'sub', sub, 'missing', up-set(out), flush=True)
        print('n', n, 'done, bad so far', bad, flush=True)
