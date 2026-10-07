# Day 227: transcription check of Jing-Liu 2104.04411 (2.32) recursion and (2.33) nested sum vs day226 green.py
import sys; sys.path.insert(0,'../day226')
import sympy as sp, functools, itertools
from green import X as Xg, parts, t, zee
def zt(r): return sp.Rational(zee(r))*sp.prod([1/(1-t**k) for k in r])
def subparts(mu):
    s=set()
    for k in range(len(mu)+1):
        for c in itertools.combinations(range(len(mu)),k):
            s.add(tuple(mu[i] for i in c))
    # multiplicity: tau ⊳ mu counts index subsets, so return list with repetition
    out=[]
    for k in range(len(mu)+1):
        for c in itertools.combinations(range(len(mu)),k): out.append(tuple(mu[i] for i in c))
    return out
@functools.lru_cache(None)
def XJL(lam,mu):   # superscript lam = HL index, subscript mu = class
    if len(lam)<=1: return sp.Integer(1)
    lam1=lam[1:]; m=sum(lam1); tot=0
    for tau in subparts(mu):
        if sum(tau)>m: continue
        for rho in parts(m-sum(tau)):
            cls=tuple(sorted(tau+rho,reverse=True))
            tot+=(-1)**len(rho)/zt(rho)*XJL(lam1,cls)
    return sp.factor(sp.simplify(tot))
if __name__=='__main__':
    bad=0;cnt=0
    for n in range(1,7):
        for lam in parts(n):
            for mu in parts(n):
                d=sp.simplify(XJL(lam,mu)-Xg(lam,mu)); cnt+=1
                if d!=0: bad+=1; print('MISMATCH',lam,mu,XJL(lam,mu),Xg(lam,mu))
    print('checked',cnt,'mismatches',bad)
