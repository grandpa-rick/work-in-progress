import sys
from fractions import Fraction as Fr
sys.path.insert(0,'../prove249c')
from eng import ell, leq, subsets
from dem import word, Di
from sympy import factorint
n=int(sys.argv[1]); P=[2,3,5,7,11,13][:n]
y=[Fr(p) for p in P]
w={m:y[n-m] for m in range(1,n+1)}   # w_m = y_{n+1-m}
def Dchar(kap,t):
    f={tuple([t]*n):1}
    for i in reversed(word(kap,n)): f=Di(f,i,n)
    s=Fr(0)
    for mu,a in f.items():
        p=Fr(a)
        for e,yy in zip(mu,y): p*=yy**e
        s+=p
    return s
def wp(M):
    p=Fr(1)
    for m in M: p*=w[m]
    return p
def mono(q):
    d={}
    for p,e in factorint(q.numerator).items(): d[p]=e
    for p,e in factorint(q.denominator).items(): d[p]=d.get(p,0)-e
    # express in w: prime P[i]=y_{i+1}=w_{n-i}
    return {('w%d'%(n-P.index(p))):e for p,e in d.items() if p in P}, (q.numerator, q.denominator) if any(p not in P for p in list(factorint(q.numerator))+list(factorint(q.denominator))) else ''
for kap in subsets(n):
    R=[sum(ell(kap,M,n,w)*wp(M)**(-2*t-2) for M in subsets(n) if leq(M,kap))/Dchar(set(kap),t) for t in (0,1)]
    print(kap,'A=',mono(R[0]),' B=',mono(R[1]/R[0]))
