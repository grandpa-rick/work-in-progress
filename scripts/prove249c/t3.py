import sympy as sp, itertools, sys
from engs import *
n=int(sys.argv[1]); c=int(sys.argv[2])
W=sp.symbols('w1:%d'%(n+1)); w={k:W[k-1] for k in range(1,n+1)}
D=sp.Integer(1)
for i in range(1,n+1):
    D*=(w[i]-1/w[i])
    for j in range(i+1,n+1): D*=(w[j]-w[i])*(w[i]*w[j]-1)
def vertices(pts):
    # vertices of Newton polytope via random linear functionals (approx: extremal points)
    import random
    V=set()
    for _ in range(3000):
        v=[random.gauss(0,1) for _ in range(n)]
        best=max(pts,key=lambda p:sum(a*b for a,b in zip(p,v)))
        V.add(best)
    return V
cases=eval(sys.argv[3])
for lam,mu,kap in cases:
    X=Cfun(list(lam),list(mu),list(kap),c,n,w)
    P=sp.cancel(sp.together(D*X))
    num,den=sp.fraction(P)
    num=sp.Poly(sp.expand(num),*W); 
    dd=sp.Poly(sp.expand(den),*W)
    print(lam,mu,kap,'den terms',len(dd.terms()), dd.terms()[:2])
    d=dict(num.terms()); V=vertices(list(d))
    print('  #terms',len(d),' vertex coeffs',sorted(set(d[v] for v in V)), '#V',len(V))
