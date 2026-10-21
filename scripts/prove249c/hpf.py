from hp import *
import sys
n=int(sys.argv[1]); a=[Fr(k*k+1,7) for k in range(1,n+1)]
U=list(range(1,n+1)); subs=[list(S) for r in range(n+1) for S in itertools.combinations(U,r)]
dim=n*(n+1)//2; cs=list(range(-2,-2+2*(dim+2),2)); c_=sp.Symbol('c')
for kap in subs:
 for lam in subs:
  for mu in subs:
    if lam>mu: continue
    if not any(leq(M,kap) and leq(lam,M) and leq(mu,M) for M in subs): continue
    pts=[(c,Plim(lam,mu,kap,c,n,a)) for c in cs]
    P=sp.expand(sp.interpolate(pts,c_))
    if P==0: continue
    print(lam,mu,kap,sp.factor(P),flush=True)
