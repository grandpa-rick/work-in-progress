import sys, itertools, sympy as sp
from pq import polys, x as X
def Hres(kap, c):
    lam=sorted(kap,reverse=True); l=len(lam)
    t=sp.symbols('t1:%d'%(l+1))
    expr=1
    for s in range(1,l+1):
        ts=t[s-1]
        rho=ts**(2*c-1-s)*sp.Mul(*[(t[r-1]-ts)/(t[r-1]*ts-1) for r in range(1,s)])
        expr*=ts*(1-rho)/((ts**2-1)*(ts-1)**lam[s-1])
    for s in range(l,0,-1):   # innermost t_l first
        expr=sp.residue(sp.together(expr),t[s-1],1)
    return sp.nsimplify(sp.simplify((-1)**(sum(lam)-l)*expr))
ok=tot=0
for kap in [(1,),(2,),(3,),(1,2),(1,3),(2,3),(1,2,3),(2,4),(1,2,4),(3,4)]:
    S,E,O,p,q=polys(kap)
    for c in [sp.Rational(1,2),1,sp.Rational(3,2),2,sp.Rational(5,2),3,sp.Rational(7,2)]:
        v=Hres(kap,c); w=E.subs(X,-2*c); tot+=1; ok+= v==w
        if v!=w: print('MISS',kap,c,v,w)
print(ok,'/',tot)
