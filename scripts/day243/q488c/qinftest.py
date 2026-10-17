"""For given triples: full-q residue on H_ab (pairs from Clio's predicted set) via the lower ideal of nu;
report q-valuation at 0 and degree at infinity vs e-1."""
import itertools, sys, sympy as sp, time
from engine import geq, subsets
from fullc import cfull, eps, q
N=int(sys.argv[1]); triples=eval(open(sys.argv[2]).read())
S=subsets(N); dim=N*(N+1)//2
def rm(M,a,b): return tuple(x for x in M if x not in (a,b))
for (l,m,nu) in triples:
    pairs=[]
    for k in S:
        if not geq(nu,k): continue
        for a,b in itertools.combinations(k,2):
            if geq(rm(k,a,b),l) and geq(rm(k,a,b),m) and (a,b) not in pairs: pairs.append((a,b))
    e=dim+sum(nu)-sum(l)-sum(m); out=[]
    for a,b in pairs:
        cv=cfull(l,m,nu,N,a,b)
        r=sp.cancel(sp.cancel(cv*eps).subs(eps,0))
        if r==0: out.append(((a,b),'none')); continue
        n_,d_=sp.fraction(r); Pn=sp.Poly(n_,q); Pd=sp.Poly(d_,q)
        v0=min(x[0] for x in Pn.monoms())-min(x[0] for x in Pd.monoms()); di=Pn.degree()-Pd.degree()
        out.append(((a,b),'val0=%d'%v0,'QINF' if di==e-1 else 'deg_inf=%d (e-1=%d)'%(di,e-1)))
    print(l,m,nu,out,flush=True)
