# Hilbert-polynomial limit P(c) = lim_{w->1} X(c)/X_triv(2n), type C
from eng import *
import sys, itertools, sympy as sp
def Plim(lam,mu,kap,c,n,a):
    vals=[]
    for e in [Fr(1,10**12),Fr(2,10**12)]:
        w={i:1+a[i-1]*e for i in range(1,n+1)}
        vals.append(Cfun(lam,mu,kap,c,n,w)/Cfun([],[],list(range(1,n+1)),2*n,n,w))
    v=2*vals[0]-vals[1]
    return v.limit_denominator(4096)
if __name__=="__main__":
    n=int(sys.argv[1]); a=[Fr(k*k+1,7) for k in range(1,n+1)]
    U=list(range(1,n+1)); subs=[list(S) for r in range(n+1) for S in itertools.combinations(U,r)]
    dim=n*(n+1)//2; cs=list(range(-2,-2+2*(dim+2),2))   # even c
    c_=sp.Symbol('c'); Umax=2*n-2
    res=[]
    for kap in subs:
     for lam in subs:
      for mu in subs:
        if lam>mu: continue
        if not any(leq(M,kap) and leq(lam,M) and leq(mu,M) for M in subs): continue
        pts=[(c,Plim(lam,mu,kap,c,n,a)) for c in cs]
        P=sp.expand(sp.interpolate(pts,c_))
        # check extra point odd
        chk=Plim(lam,mu,kap,2*n+3,n,a); ok=(sp.Rational(chk.numerator,chk.denominator)==P.subs(c_,2*n+3))
        if P==0: res.append((lam,mu,kap,'P==0')); print(lam,mu,kap,'P identically 0'); continue
        rts=[r for r in sp.Poly(P,c_).nroots()] 
        mr=max(sp.re(r) for r in rts) if rts else None
        realr=sorted(round(float(sp.re(r)),3) for r in rts if abs(sp.im(r))<1e-9)
        print(lam,mu,kap,'deg',sp.degree(P,c_),'oddchk',ok,'maxRe root',round(float(mr),3) if mr is not None else None,'real roots',realr, 'P(U+1)',P.subs(c_,Umax+1))
