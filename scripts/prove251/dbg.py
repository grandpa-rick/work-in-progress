import sympy as sp, itertools
from fractions import Fraction as Fr
import pf
from dectest import JS
lam=[3,2,1]; c=sp.Rational(5,2); N=5
def Jz(S):
    S=list(S); pairs=[(r,s) for i,r in enumerate(S) for s in S[i+1:]]; maxn=sum(lam)
    tot=Fr(0)
    def rec(pi,shift,coef):
        nonlocal tot
        if pi==len(pairs):
            v=coef
            for s in S: v*=pf.coeff_f(N,s,lam[s-1],lam[s-1]-shift[s])*Fr(1,2**(lam[s-1]+1))
            tot+=v; return
        r,s=pairs[pi]
        for n in range(maxn+1):
            sh=dict(shift); sh[r]-=n; sh[s]+=n
            rec(pi+1,sh,coef*(1 if n==0 else 2*(-1)**n))
    rec(0,{s:0 for s in S},Fr(1)); return tot
for q in range(1,4):
    for S in itertools.combinations([1,2,3],q):
        print(S, JS(lam,c,S), Jz(S))
