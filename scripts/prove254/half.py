# half-integer values s_{a+1/2}(kappa) via (*) from anchor s_{1/2}=(-1)^{|k|}2^{-k1}, |c|=V(theta) (Q30)
import itertools, sys
from fractions import Fraction as Fr
from functools import lru_cache
def boxes(k): return frozenset((i,i+j) for i,li in enumerate(k) for j in range(li))
def sp(n): return [tuple(sorted(c,reverse=True)) for L in range(n+1) for c in itertools.combinations(range(1,n+1),L)]
def le(a,b): return a[0]<=b[0] and a[1]<=b[1]
def V(th):
    th=list(th)
    for x in th:
        for y in th:
            if x!=y and le(x,y):
                if x[0]==x[1] or y[0]==y[1] and False: pass
    # closed form: no 3-chain, no diagonal non-minimal; 2^{isolated off-diag}
    S=set(th)
    for x in th:
        for y in th:
            if x!=y and le(x,y):
                if y[0]==y[1]: return 0
                for z in th:
                    if z!=y and le(y,z): return 0
    iso=sum(1 for x in th if x[0]!=x[1] and not any(y!=x and (le(x,y) or le(y,x)) for y in th))
    return 2**iso
N=int(sys.argv[1]) if len(sys.argv)>1 else 5
KS=sorted(sp(N),key=sum)
SUP={k:{r:V(boxes(k)-boxes(r)) for r in KS if boxes(r)<boxes(k)} for k in KS}
SUP={k:{r:v for r,v in d.items() if v} for k,d in SUP.items()}
@lru_cache(None)
def s(k,x2):  # x = x2/2
    if not k: return Fr(1)
    if x2==1: return Fr((-1)**sum(k),2**k[0])
    if x2%2==0 and x2==2: return Fr(0)
    y=x2-2
    return s(k,y)+sum(v*s(r,y) for r,v in SUP[k].items())
if __name__=='__main__':
  for k in KS:
    if k: print(k, [s(k,2*a+1)*2**k[0] for a in range(0,7)], 'int', [s(k,2*j) for j in range(1,6)])
