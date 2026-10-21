# lowest epsilon-order term P_d(c;a) of X(c)/X_triv(2n) at w=1+a eps, type C
from eng import *
import sys, itertools, numpy as np
def val(lam,mu,kap,c,n,a,e):
    w={i:1+a[i-1]*e for i in range(1,n+1)}
    return Cfun(lam,mu,kap,c,n,w)/Cfun([],[],list(range(1,n+1)),2*n,n,w)
def lead(lam,mu,kap,c,n,a,d=None):
    e1=Fr(1,10**15); 
    v1=val(lam,mu,kap,c,n,a,e1); v2=val(lam,mu,kap,c,n,a,2*e1); v4=val(lam,mu,kap,c,n,a,4*e1)
    if d is None:
        if v1==0: return None,0
        import math
        d=round(math.log2(abs(float(v2/v1))))
    # f = e^d (P + Q e + ...): Richardson
    g1=v1/e1**d; g2=v2/(2*e1)**d; g4=v4/(4*e1)**d
    # second order extrapolation
    P=(8*g1-6*g2+g4)/3
    return d,float(P)
if __name__=="__main__":
    n=int(sys.argv[1]); seed=int(sys.argv[2])
    import random; random.seed(seed)
    a=[Fr(random.randint(1,20),random.randint(1,20))*random.choice([1,-1]) for k in range(n)]
    U=list(range(1,n+1)); subs=[list(S) for r in range(n+1) for S in itertools.combinations(U,r)]
    dim=n*(n+1)//2; Umax=2*n-2
    cs=list(range(Umax+1-2*(dim+3),Umax+2+2*(dim+3)))  # wide integer range
    worst=-99; out=[]
    for kap in subs:
     for lam in subs:
      for mu in subs:
        if lam>mu: continue
        if not any(leq(M,kap) and leq(lam,M) and leq(mu,M) for M in subs): continue
        # find d at a large c
        d,_=lead(lam,mu,kap,Umax+7,n,a)
        vals=[lead(lam,mu,kap,c,n,a,d)[1] for c in cs]
        co=np.polyfit(cs,vals,dim+2)
        # trim tiny leading coeffs
        sc=max(abs(v) for v in vals)
        co=np.trim_zeros(np.where(abs(co)*max(abs(np.array(cs)))**np.arange(len(co))[::-1]<1e-7*sc,0,co),'f')
        r=np.roots(co) if len(co)>1 else []
        mr=max([x.real for x in r],default=None)
        realr=sorted(round(x.real,3) for x in r if abs(x.imag)<1e-6)
        intabove=[c for c in range(Umax+1,Umax+40) if abs(np.polyval(co,c))<1e-6*sc]
        print(lam,mu,kap,'d',d,'deg',len(co)-1,'maxRe',None if mr is None else round(mr,3),'real',realr,'intzeros>U',intabove,flush=True)
