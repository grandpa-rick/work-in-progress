"""Unpadded evaluation: A_{M,I} = WC~_{N,s,r}(w_M|w|I), s=|M| (Motegi stability eq:typeC-stability-s-ge-r).
Returns A0=A|_{q=0}, A1=dA/dq|_{q=0}. Type C only."""
import itertools, sympy as sp
from engine import summand, denom_pairs, geq
q=sp.Symbol('q')
def A01(N,M,I):
    s,r=len(M),len(I)
    if not geq(M,I): return sp.Integer(0),sp.Integer(0)
    w=sp.symbols('w1:%d'%(N+1)); n=s
    u=[w[m-1] for m in M]
    tot=0
    for perm in itertools.permutations(range(n)):
        us=[u[p] for p in perm]
        t=summand(us,w,I,N,n,1,0,q)
        d=sp.prod([us[a]-us[b] for (a,b) in denom_pairs(n,r)])
        tot+=t/d
    val=tot/sp.factorial(n-r)*sp.prod([x**(-(N+n)) for x in u])
    E=(1-q)**(-(n-r))*sp.prod([1-q*a/b for a in u for b in u])*sp.prod([1-q*a*b for a in u for b in w])
    E/=sp.prod([1-q*u[j]*u[k] for j in range(n) for k in range(j+1,n)])
    E/=sp.prod([1-q*w[j]*w[k] for j in range(N) for k in range(j,N)])
    val=val/E
    a0=sp.factor(sp.cancel(val.subs(q,0)))
    a1=sp.factor(sp.cancel(sp.diff(val,q).subs(q,0)))
    return a0,a1
