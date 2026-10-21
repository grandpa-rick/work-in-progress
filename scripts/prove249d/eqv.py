# Equivariant test: ell_kappa(w^c) vs Demazure character D_w(e^{t omega}), c=-2t-2, several matchings
import sys, itertools, random
from fractions import Fraction as Fr
sys.path.insert(0,'../prove249c')
from eng import ell, leq, subsets
from dem import word, Di
n=int(sys.argv[1]); random.seed(int(sys.argv[2]) if len(sys.argv)>2 else 1)
y=[Fr(random.randint(2,30),random.randint(2,30)) for _ in range(n)]
def Dchar(kap,t):
    f={tuple([t]*n):1}
    for i in reversed(word(kap,n)): f=Di(f,i,n)
    s=Fr(0)
    for mu,a in f.items():
        p=Fr(a)
        for e,yy in zip(mu,y): p*=yy**e
        s+=p
    return s
def wp(w,M):
    p=Fr(1)
    for m in M: p*=w[m]
    return p
maps={'rev':lambda m:y[n-m],'revinv':lambda m:1/y[n-m],'id':lambda m:y[m-1],'idinv':lambda m:1/y[m-1]}
for kap in subsets(n):
    if not kap: continue
    res=[]
    for nm,f in maps.items():
        w={m:f(m) for m in range(1,n+1)}
        R=[]
        for t in range(0,5):
            c=-2*t-2
            L=sum(ell(kap,M,n,w)*wp(w,M)**c for M in subsets(n) if leq(M,kap))
            R.append(L/Dchar(set(kap),t))
        q=[R[i+1]/R[i] for i in range(len(R)-1)]
        res.append((nm, len(set(q))==1))
    print(kap,res)
