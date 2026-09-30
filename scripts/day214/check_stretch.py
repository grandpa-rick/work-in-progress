"""Day 214 stretch checks.
(1) identity sum_{|B|=r} prod_{i in B, j notin B} x_i/(x_i-x_j) = 1 (N<=5) and the t-version = [N choose r]_t (N<=4).
(2) conjugate peel formula: peel_k(kappa)'_c = kappa'_c - min(k,kappa'_c) + min(k,kappa'_{c+1}).
(3) at t=0: c_{lam mu}(s,0) = s^{n(mu)} (1 + O(s)) for all mu ⊵ lam, |lam| <= N (exact interpolation in s with check point);
    at generic t: val_s c_{lam mu} = n(mu)."""
import sympy as sp, itertools, sys
from fractions import Fraction as Fr
from ds_pointeval import parts, dominates, nstat, coeffs, solve
import random
ok = True
t = sp.symbols('t')
for N in range(1, 6):
    xs = sp.symbols(f'y1:{N+1}')
    for r in range(0, N+1):
        S0 = sum(sp.prod([xs[i]/(xs[i]-xs[j]) for i in B for j in range(N) if j not in B]) for B in itertools.combinations(range(N), r))
        ok &= sp.cancel(sp.together(S0 - 1)) == 0
        if N <= 4:
            St = sum(sp.prod([(xs[i]-t*xs[j])/(xs[i]-xs[j]) for i in B for j in range(N) if j not in B]) for B in itertools.combinations(range(N), r))
            tb = sp.Integer(1)*sp.prod([1-t**(N-i) for i in range(r)], sp.Integer(1))/sp.prod([1-t**(i+1) for i in range(r)], sp.Integer(1))
            ok &= sp.cancel(sp.together(St - tb)) == 0
print('(1) level-set kernel identities:', ok, flush=True)
def conj(k):
    return tuple(sum(1 for p in k if p > i) for i in range(k[0])) if k else ()
def peel(kap, k): return tuple(sorted((x for x in [kap[i]-1 if i < k else kap[i] for i in range(len(kap))] if x > 0), reverse=True))
ok2 = True
for n in range(1, 13):
    for kap in parts(n):
        kc = list(conj(kap)) + [0]
        for k in range(1, len(kap)+1):
            pred = [kc[c] - min(k, kc[c]) + min(k, kc[c+1]) for c in range(len(kc)-1)]
            pred = tuple(x for x in pred if x > 0)
            ok2 &= pred == conj(peel(kap, k))
print('(2) conjugate peel formula:', ok2, flush=True)
ok &= ok2
N = int(sys.argv[1]); rng = random.Random(7)
for T in (Fr(0), Fr(-5, 11)):
    for n in range(2, N+1):
        for lam in parts(n):
            D = nstat(lam) + 1
            sv = [Fr(i+2, 3) for i in range(D+1)]
            data = [coeffs(lam, x, T, rng) for x in sv]
            chk = coeffs(lam, Fr(17, 5), T, rng)
            for mu in parts(n):
                if not dominates(mu, lam): continue
                c = solve([[x**j for j in range(D+1)] for x in sv], [d.get(mu, Fr(0)) for d in data])
                assert sum(cj*Fr(17, 5)**j for j, cj in enumerate(c)) == chk.get(mu, 0), ('degree bound', lam, mu)
                v = next(i for i, a in enumerate(c) if a != 0)
                good = v == nstat(mu) and (T != 0 or c[v] == 1)
                if not good: print('FAIL', T, lam, mu, v, c[v])
                ok &= good
        print('(3) t =', T, 'n =', n, 'done', flush=True)
print('ALL OK' if ok else 'FAIL')
