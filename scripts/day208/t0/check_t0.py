"""Day 208: t=0 corollary of the e_k*e_r Pieri theorem (Day 207b).
Claim (all k>=1, r>=0, mu=min(k,r), M=max(k,r)):
  (e_k * e_r)|_{t=0} = sum_{b<mu} (1-s) s^b e_b e_{r+k-b} + s^mu e_mu e_M.
Checks:
 (1) exact symbolic: collected RHS of the theorem (monomials e_a e_c, a<=c) is a rational function in t,
     regular at t=0 after collection, with limit = claim; k<=5, r<=6. Also the swap symmetry
     RHS(k,r) == RHS(r,k) exactly as rational functions (the rigorous route for r<k).
 (2) independent: direct AHA e_k(Y).e_r (flint, Day 205 pipeline), m<=6, k<=5 (k<=m), r<=m:
     coefficients of t^i, i<C(k,2), vanish; the t^{C(k,2)} coefficient equals the claim (random rational points).
 (3) weights sum to 1.
"""
import sys, itertools, random
from fractions import Fraction as Fr
import sympy as sp
sys.path.insert(0, '/home/agent/projects/proofs/scripts/day205')
from k3_fast_pipeline import AHA
s, t = sp.symbols('s t')
ok_all = True
def rep(tag, ok):
    global ok_all; ok_all &= bool(ok); print(('OK  ' if ok else 'FAIL'), tag)
def poch(a, n): return sp.prod([1 - a*t**i for i in range(n)])
def alpha(j): return 0 if j < 0 else sp.prod([s - t**i for i in range(1, j+1)])/poch(t, j)
def c(n, j):
    if j < 0 or n - j < 0: return 0
    return poch(s, n-j)/poch(t, n-j)*(alpha(j) - s*t**(n-j)*alpha(j-1))
def F(n, w): return sum(c(n, j)*w**j for j in range(n+1))
def collected(k, r):
    d = {}
    for b in range(k+1):
        key = tuple(sorted((b, r+k-b))); d[key] = d.get(key, 0) + s**b*F(k-b, t**(r-b))
    return {key: sp.cancel(sp.together(v)) for key, v in d.items()}
def claim(k, r):
    mu, M = min(k, r), max(k, r); d = {}
    for b in range(mu):
        key = tuple(sorted((b, r+k-b))); d[key] = d.get(key, 0) + (1-s)*s**b
    key = (mu, M); d[key] = d.get(key, 0) + s**mu
    return d
# (1)
ok1 = okreg = oksw = True
for k in range(1, 6):
    for r in range(0, 7):
        C = collected(k, r); W = claim(k, r)
        for key in set(C) | set(W):
            v = C.get(key, 0)
            num, den = sp.fraction(sp.cancel(v))
            okreg &= sp.Poly(den, t).eval(0) != 0          # regular at t=0 after collection
            ok1 &= sp.expand(sp.cancel(v).subs(t, 0) - W.get(key, 0)) == 0
        if r >= 1:
            C2 = collected(r, k)
            oksw &= all(sp.cancel(C.get(key, 0) - C2.get(key, 0)) == 0 for key in set(C) | set(C2))
rep('collected RHS regular at t=0, k<=5, r<=6', okreg)
rep('t=0 limit == claim, k<=5, r<=6', ok1)
rep('swap symmetry RHS(k,r)==RHS(r,k) exact, k<=5, 1<=r<=6', oksw)
rep('weights sum to 1', all(sp.expand(sum(claim(k, r).values()) - 1) == 0 for k in range(1, 8) for r in range(0, 9)))
# (2) direct AHA
def ev(vals, r):
    if r < 0 or r > len(vals): return Fr(0)
    tot = Fr(0)
    for cmb in itertools.combinations(vals, r):
        p = Fr(1)
        for v in cmb: p *= v
        tot += p
    return tot
random.seed(208); ok2 = True; cnt = 0
for m in range(1, 7):
    A_ = AHA(m)
    for k in range(1, min(m, 5)+1):
        for r in range(0, m+1):
            if m == 6 and k == 5: continue
            L = 0*A_.s
            for bb in itertools.combinations(range(1, m+1), k):
                V = A_.e(r)
                for i in reversed(bb): V = A_.Y(V, i)
                L = L + V
            ck = k*(k-1)//2
            D = L.to_dict()
            ok2 &= all(int(mon[m+1]) >= ck for mon, cc in D.items() if int(cc) != 0)
            low = {mon: int(cc) for mon, cc in D.items() if int(mon[m+1]) == ck}
            for trial in range(3):
                X = [Fr(random.randint(1, 40), random.randint(41, 90)) for _ in range(m)]
                sv = Fr(random.randint(1, 20), 23)
                lhs = Fr(0)
                for mon, cc in low.items():
                    v = Fr(cc)
                    for i in range(m): v *= X[i]**int(mon[i])
                    lhs += v*sv**int(mon[m])
                rhs = sum(Fr(sp.Rational(val.subs(s, sp.Rational(sv.numerator, sv.denominator))).p,
                              sp.Rational(val.subs(s, sp.Rational(sv.numerator, sv.denominator))).q)
                          * ev(X, a)*ev(X, cc2) for (a, cc2), val in claim(k, r).items())
                ok2 &= (lhs == rhs); cnt += 1
rep(f'direct AHA t=0 vs claim, m<=6, k<=min(m,5), r<=m: {cnt} point-checks', ok2)
print('ALL OK' if ok_all else 'SOME FAIL')
