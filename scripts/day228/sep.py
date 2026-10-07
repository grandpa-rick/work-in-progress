"""Day 228 Item 1: separator in ONE normalisation.
Gauge group on full-merge leads L(lam): L -> L * prod f(lam_i) * g(|lam|) * a^{l(lam)}  (f,g,a arbitrary).
I  = L(111)/(L(21)L(11))      : invariant under f, a^{l-1}; NOT under g (changes by 1/g(2)).
J  = L(1111)L(22)/L(211)^2    : invariant under f, g, a^l.
Objects: Rick (Thm G formula; checked vs raw star pickles earlier), HT phi-reading <kappa_HT(lam),e_n>,
Dolega (b) column-cumulant [g_n] lead (dict_test.py)."""
import sys; sys.path.insert(0,'/home/agent/projects/scripts/browse_dolega'); sys.path.insert(0,'/home/agent/projects/scripts/day216b')
import sympy as sp
from math import factorial
from dict_test import setparts, conn_graph_sum, RickLead, lowest
from symf import *
T=sp.Symbol('t')
def HTphi(lam):  # <kappa_HT(lam), e_n> with q:=t ; phi(h~_m)=t^{C(m,2)} (proved: e_m(1,q,q^2,..)(q;q)_m)
    r=len(lam); s=0
    for pi in setparts(range(r)):
        k=len(pi); s+=(-1)**(k-1)*factorial(k-1)*sp.prod([T**sp.binomial(sum(lam[i] for i in B),2) for B in pi])
    return sp.factor(s/(T-1)**(r-1))
def Dolb(lam):
    n=sum(lam); lp=transpose(tuple(sorted(lam,reverse=True)))
    H=Htilde(lp); gb={m:prod([Htilde((1,)*i) for i in m]) for m in parts(n)}
    d=to_basis(H,gb); v,c=lowest(d[(n,)]); assert v==len(lam)-1,(lam,v); return sp.factor(c.subs(t,T))
objs={'Rick':lambda l:RickLead(l).subs(t,T),'HTphi':HTphi,'Dol(b)':Dolb}
lams=[(1,1),(2,1),(1,1,1),(2,2),(2,1,1),(1,1,1,1)]
V={}
for name,F in objs.items():
    V[name]={l:sp.factor(F(l)) for l in lams}
    print(name,V[name])
for name in objs:
    L=V[name]
    I=sp.factor(L[(1,1,1)]/(L[(2,1)]*L[(1,1)])); J=sp.factor(L[(1,1,1,1)]*L[(2,2)]/L[(2,1,1)]**2)
    print(f"{name:7s} I={I}   J={J}   J(0)={sp.nsimplify(J.subs(T,0))}  J(1)={sp.limit(J,T,1)}")
