# configuration-level leading terms for type C problem
import itertools, random
from collections import defaultdict
def leq(A,B):
    A=sorted(A,reverse=True); B=sorted(B,reverse=True)
    return len(A)<=len(B) and all(A[i]<=B[i] for i in range(len(A)))
class Term:
    __slots__=('s','e')
    def __init__(s,n): s.s=1; s.e=[0]*(n+1)
def addm(T,m,k=1):
    for i,v in m.items(): T.e[i]+=k*v
def fac(T,m,p,xi):
    # multiply by (1 - w^m)^p, leading term under xi
    d=sum(xi[i]*v for i,v in m.items())
    if d>0:
        T.s*=(-1); addm(T,m,p)   # (−m)^p : sign (−1)^p = −1 for p=±1
def T0terms(kap,M,th,xi,n):
    """list of (sign, exp) leading terms of Psi_{kap,th}(M) configurations"""
    out=[]
    r=len(kap)
    for h in itertools.permutations(M,r):
        if any(h[j]<kap[j] for j in range(r)): continue
        T=Term(n); tails=[m for m in M if m not in h]
        for j in range(r):
            x=h[j]
            fac(T,{x:-2},1,xi)
            for k in range(1,kap[j]): fac(T,{k:1,x:-1},1,xi)
            T.e[x]+=th[j]
        for j in range(r):
            for k in range(j+1,r):
                fac(T,{h[j]:-1,h[k]:-1},1,xi); fac(T,{h[j]:1,h[k]:-1},-1,xi)
        for t in tails:
            for j in range(r):
                fac(T,{t:-1,h[j]:-1},1,xi); fac(T,{t:1,h[j]:-1},-1,xi)
        out.append((T.s,tuple(T.e)))
    return out
def Fterms(kap,M,xi,n):
    r=len(kap); res=[]
    for s,e in T0terms(kap,M,[kap[j]-2*(j+1) for j in range(r)],xi,n):
        e=list(e)
        for m in M: e[m]+=2*r
        res.append((s,tuple(e)))
    return res
def single(ts):
    assert len(ts)==1; return ts[0]
def ctrm(I,n):
    s=1; e=[0]*(n+1)
    for i in I:
        s*=(-1)**i
        for k in range(1,i+1): e[k]+=1
    return s,e
def Cconfigs(lam,mu,kap,c,n,xi):
    U=list(range(1,n+1)); kc=[i for i in U if i not in kap]
    out=[]
    for r in range(n+1):
      for M in itertools.combinations(U,r):
        M=list(M)
        if not (leq(M,kap) and leq(lam,M) and leq(mu,M)): continue
        Mc=[i for i in U if i not in M]
        s0,e0=ctrm(M,n); s1,e1=ctrm(kap,n)
        sF,eF=single(Fterms(M,M,xi,n)); sP,eP=single(T0terms(Mc,Mc,[0]*len(Mc),xi,n))
        base_s=s0*s1*sF*sP
        base=[e0[i]-e1[i]-eF[i]-eP[i]+(c if i in M else 0) for i in range(n+1)]
        for sd,ed in T0terms(kc,Mc,[0]*len(kc),xi,n):
          for sl,el in Fterms(lam,M,xi,n):
            for sm,em in Fterms(mu,M,xi,n):
              out.append((base_s*sd*sl*sm, tuple(base[i]+ed[i]+el[i]+em[i] for i in range(n+1)), tuple(M)))
    return out
def top(cfgs,xi):
    g=defaultdict(int); cnt=defaultdict(int); Ls=defaultdict(set)
    for s,e,M in cfgs: g[e]+=s; cnt[e]+=1; Ls[e].add(M)
    es=sorted(g,key=lambda e:-sum(xi[i]*e[i] for i in range(len(e))))
    return [(e,g[e],cnt[e],Ls[e]) for e in es]
