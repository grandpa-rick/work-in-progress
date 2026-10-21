# Expand F_kappa(M) * w_M^a in basis G_nu(M) (G = Psi_{nu,0} or its iota), report coefficients as monomials at prime point
import sys
from fractions import Fraction as Fr
sys.path.insert(0,'../prove249c')
from eng import F, Psi, leq, subsets
from sympy import factorint, Matrix, Rational
n=int(sys.argv[1]); a=int(sys.argv[2]); mode=sys.argv[3]
P=[2,3,5,7,11,13][:n]; w={i:Fr(P[i-1]) for i in range(1,n+1)}; wi={i:1/w[i] for i in w}
S=[tuple(s) for s in subsets(n)]
def wp(M,ww):
    p=Fr(1)
    for m in M: p*=ww[m]
    return p
def G(v,M):
    if mode=='psi': return Psi(list(v),[0]*len(v),list(M),w)
    else: return Psi(list(v),[0]*len(v),list(M),wi)
def mono(q):
    if q==0: return 0
    d={}
    s='-' if q<0 else ''
    for p,e in factorint(abs(q.numerator)).items(): d[p]=e
    for p,e in factorint(q.denominator).items(): d[p]=d.get(p,0)-e
    bad=[p for p in d if p not in P]
    return s+' '.join('w%d^%d'%(P.index(p)+1,e) for p,e in sorted(d.items()) if p in P)+(' *junk%s'%bad if bad else '')
Gm=Matrix([[Rational(*G(v,M).as_integer_ratio()) for v in S] for M in S])
Gi=Gm.inv()
for k in S:
    f=Matrix([Rational(*(F(list(k),list(M),w)*wp(M,w)**a).as_integer_ratio()) for M in S])
    co=Gi*f
    print(k,{S[i]:mono(Fr(int(co[i].p),int(co[i].q))) for i in range(len(S)) if co[i]!=0})
