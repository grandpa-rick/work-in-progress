import sys; sys.path.insert(0,'../prove247')
import mp
mp.K=14
from mp import S,Sys,subsets,geq,P,inv
def ell(s,nu,L): return s.coord({L:S.const(1)},nu)
for n in (4,5):
  for nu in subsets(n):
    if n not in nu or len(nu)==n: continue
    kap=[k for k in range(1,n+1) if k not in nu]; m=max(kap)
    if len(kap)==1 or kap[-2]<m-1: continue
    s=Sys(n,7); s.w[n-1]=S(1,[1]+[0]*(mp.K-1)); t=Sys(n-1,7)
    for i in range(n-1): t.w[i]=s.w[i]
    Ls=[L for L in subsets(n-1) if geq(nu,L)]
    lam={}
    for L in Ls:
        v=ell(s,nu,L).norm(); lam[L]=v.coef(n-m-1)
    hits=[]
    for nup in subsets(n-1):
      for e in range(-4,5):
        rats=set(); bad=False
        for L in subsets(n-1):
            ref=ell(t,nup,L) if geq(nup,L) else S.const(0)
            ref=ref*t.wpow(L,e) if not ref.iszero() else ref
            lv=lam.get(L,0)
            if ref.iszero():
                if lv!=0: bad=True;break
            else: rats.add(lv*inv(ref.c[0])%P)
        if not bad and len(rats)==1 and 0 not in rats: hits.append((nup,e))
    print(n,nu,'kap',kap,'hits',hits)
