import sympy as sp, itertools, sys, random
from engs import *
n=int(sys.argv[1]); c=int(sys.argv[2])
W=sp.symbols('w1:%d'%(n+1)); w={k:W[k-1] for k in range(1,n+1)}
D=sp.Integer(1)
for i in range(1,n+1):
    D*=(w[i]-1/w[i])
    for j in range(i+1,n+1): D*=(w[j]-w[i])*(w[i]*w[j]-1)
def lp(expr):
    P=sp.cancel(sp.together(expr)); num,den=sp.fraction(P)
    dd=sp.Poly(sp.expand(den),*W).terms(); assert len(dd)==1
    de,dc=dd[0]
    return {tuple(a-b for a,b in zip(e,de)):co/dc for e,co in sp.Poly(sp.expand(num),*W).terms()}
def verts(d):
    pts=list(d); V=set()
    for _ in range(4000):
        v=[random.gauss(0,1) for _ in range(n)]
        V.add(max(pts,key=lambda p:sum(a*b for a,b in zip(p,v))))
    return V
cases=eval(sys.argv[3])
for lam,mu,kap in cases:
    d=lp(D*Cfun(list(lam),list(mu),list(kap),c,n,w)); V=verts(d)
    # shift: center
    cen=[sp.Rational(sum(v[i] for v in V),len(V)) for i in range(n)]
    rel=sorted(set(tuple(sorted(abs(v[i]-cen[i]) for i in range(n))) for v in V))
    print(lam,mu,kap,'#V',len(V),'center',cen,'|v-cen| sorted types',rel, 'coeffs',sorted(set(d[v] for v in V)))
