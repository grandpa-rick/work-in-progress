# left W-action on g_lam(L)=w_L^2 F_lam(L): (s.g)(L) = g|_{w->s w}(s L); expand in basis g_nu
import sys, random
from fractions import Fraction as Fr
sys.path.insert(0,'../prove249c')
from eng import F, leq, subsets
from sympy import Matrix, Rational, factorint
n=int(sys.argv[1]); P=[int(x) for x in sys.argv[2].split(",")][:n]
w={i:Fr(P[i-1]) for i in range(1,n+1)}
S=[tuple(s) for s in subsets(n)]
def sw(i,w):
    w=dict(w)
    if i<n: m=n-i; w[m],w[m+1]=w[m+1],w[m]
    else: w[1]=1/w[1]
    return w
def sL(i,L):
    L=set(L)
    if i<n:
        m=n-i; L={ (m+1 if x==m else m if x==m+1 else x) for x in L}
    else: L=L^{1}
    return tuple(sorted(L))
def g(lam,L,ww):
    p=Fr(1)
    for m in L: p*=ww[m]**2
    return p*F(list(lam),list(L),ww)
def mono(q):
    if q==0: return '0'
    d={}
    for p,e in factorint(abs(q.numerator)).items(): d[p]=e
    for p,e in factorint(q.denominator).items(): d[p]=d.get(p,0)-e
    bad=[p for p in d if p not in P]
    return ('-' if q<0 else '')+' '.join('w%d^%d'%(P.index(p)+1,e) for p,e in sorted(d.items()) if p in P)+(' [nonmono %s]'%q if bad else '')
G=Matrix([[Rational(*g(v,M,w).as_integer_ratio()) for v in S] for M in S]); Gi=G.inv()
for i in range(1,n+1):
    for lam in S:
        sw_=sw(i,w)
        f=Matrix([Rational(*g(lam,sL(i,M),sw_).as_integer_ratio()) for M in S])
        co=Gi*f
        print('s%d'%i,lam,'->',{S[k]:mono(Fr(int(co[k].p),int(co[k].q))) for k in range(len(S)) if co[k]!=0})
