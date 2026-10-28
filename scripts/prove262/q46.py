# PROVE 262 Q46: rank of {j -> spec ch_rho(j) : rho ⊆ kappa} under x_i -> q^{a_i}
# Rigorous rank: spec ch_rho(j) = sum_s A_{rho,s}(q) t_s^j with t_s=q^s, s in S={(eps,a): eps sign vector}.
# Evaluating at |S| consecutive j is injective on span{t_s^j} (Vandermonde), so rank of rows j=0..|S|-1 = true rank.
import sys, itertools
import sympy as sp
from sympy.polys.matrices import DomainMatrix
sys.path.insert(0,'prove249d'); from dem import Di, word
q=sp.Symbol('q')
def strict(m): return [tuple(sorted(c,reverse=True)) for L in range(m+1) for c in itertools.combinations(range(1,m+1),L)]
def sub(r,k): return len(r)<=len(k) and all(r[i]<=k[i] for i in range(len(r)))
_c={}
def chw(k,j,m):
    key=(k,j,m)
    if key not in _c:
        f={tuple([j]*m):1}
        for i in reversed(word(set(k),m)): f=Di(f,i,m)
        _c[key]=f
    return _c[key]
def spec(f,a,qv=q):
    return sum(c*qv**sum(e[t]*a[t] for t in range(len(a))) for e,c in f.items())
from fractions import Fraction
def fspec(f,a,qv):
    v=sum(c*Fraction(qv)**sum(e[t]*a[t] for t in range(len(a))) for e,c in f.items()); return (v.numerator,v.denominator)
def rank_for(k,a,m,qv=None):
    R=[r for r in strict(m) if sub(r,k)]
    S=sorted({sum(e*x for e,x in zip(eps,a)) for eps in itertools.product((1,-1),repeat=m)})
    rows=list(range(len(S)))
    if qv is None:
        F=sp.QQ.frac_field(q)
        M=DomainMatrix([[F.convert(sp.expand(spec(chw(r,j,m),a))) for r in R] for j in rows],(len(rows),len(R)),F)
    else:
        M=DomainMatrix([[sp.QQ(*fspec(chw(r,j,m),a,qv)) for r in R] for j in rows],(len(rows),len(R)),sp.QQ)
    return M.rank(), len(R), R
if __name__=='__main__':
    m=int(sys.argv[1]); mode=sys.argv[2] if len(sys.argv)>2 else 'rho'
    a={'rho':list(range(m,0,-1)),'one':[0]*m}.get(mode) or [int(x) for x in mode.split(',')]
    for k in strict(m):
        if mode=='one':
            # x=1: specialised values are just Euler chars; rank over Q of polys in j: use rows j=0..2|R|
            R=[r for r in strict(m) if sub(r,k)]
            M=sp.Matrix([[sum(chw(r,j,m).values()) for r in R] for j in range(m*m+3)]); print(k,len(R),M.rank()); continue
        rk,n,R=rank_for(k,a,m,qv=None if '--sym' in sys.argv else 7)
        print(k,'|R|=',n,'rank=',rk,'' if rk==n else '  <<< DEFICIENT')
