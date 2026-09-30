# check remark: d_{lam mu}(t=1) = [s^{n(mu)}] c_{lam mu}(s,1) equals M_{lam, mu'} (# 0-1 matrices, row sums lam, col sums mu')
from fractions import Fraction as Fr
from ds_pointeval import parts, dominates, nstat, coeffs, solve
import random, itertools
def conj(k): return tuple(sum(1 for p in k if p > i) for i in range(k[0])) if k else ()
def M01(rows, cols):
    if not rows: return int(all(c == 0 for c in cols))
    r, tot = rows[0], 0
    for S in itertools.combinations(range(len(cols)), r):
        if all(cols[j] > 0 for j in S): tot += M01(rows[1:], [c-1 if j in S else c for j, c in enumerate(cols)])
    return tot
rng = random.Random(3); ok = True
for n in range(2, 6):
    for lam in parts(n):
        D = nstat(lam) + 1; sv = [Fr(i+2, 3) for i in range(D+1)]
        data = [coeffs(lam, x, Fr(1), rng) for x in sv]
        for mu in parts(n):
            if not dominates(mu, lam): continue
            c = solve([[x**j for j in range(D+1)] for x in sv], [d.get(mu, Fr(0)) for d in data])
            v = next(i for i, a in enumerate(c) if a != 0)
            good = v == nstat(mu) and c[v] == M01(list(lam), list(conj(mu)))
            ok &= good
            if not good: print('FAIL', lam, mu, v, c[v], M01(list(lam), list(conj(mu))))
print('t=1 remark:', 'ALL OK' if ok else 'FAIL')
