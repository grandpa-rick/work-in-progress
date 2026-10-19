import sys; sys.path.insert(0,'../prove247')
import mp
mp.K=14
from mp import S,Sys,subsets,geq,P
import random
def lead(x):
    x=x.norm(); return x.v, x.c[0]
def ell(s,nu,L):
    return s.coord({L:S.const(1)},nu)
for n in (4,5):
  for nu in subsets(n):
    if n not in nu or len(nu)==n: continue
    kap=[k for k in range(1,n+1) if k not in nu]; m=max(kap)
    clean = len(kap)==1 or kap[-2]<m-1
    s=Sys(n,7); s.w[n-1]=S(1,[1]+[0]*(mp.K-1))   # u = eps
    t=Sys(n-1,7)
    for i in range(n-1): t.w[i]=s.w[i]
    nut=tuple(sorted([x for x in nu if x!=n]+[m]))
    rats=set(); ords=set()
    for L in subsets(n-1):
        if not geq(nu,L): continue
        v=ell(s,nu,L)
        if v.iszero(): ords.add('Z'); continue
        o,c0=lead(v); ords.add(o)
        r=ell(t,nut,L) if geq(nut,L) else S.const(0)
        if o==n-m-1: rats.add(c0*mp.inv(r.c[0])%P if not r.iszero() else 'refzero')
    print(n,nu,'kap',kap,'clean' if clean else 'DIRTY','pred',n-m-1,'ords',sorted(ords,key=str),'nrat',len(rats))
