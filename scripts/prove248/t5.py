from psi import *
from t2 import lead0
import sys
# Lemma P: Psi_{(s..s+r-1),t}(M), u=w_n top, heads lead: (-1)^{s+r-1} W_{<s} prod_{M-} w^-2 u^{t-s-r} Psi_{(s..s+r-2),t+1}(M-)
for s,r,t,M in [(2,3,0,(2,3,4,5,6)),(1,2,0,(1,2,3)),(2,2,1,(2,3,4,5)),(3,2,0,(3,4,5)),(2,2,-1,(2,3,5)),(1,3,0,(1,2,3,4))]:
    n=M[-1]; u=w(n); Mm=M[:-1]
    kap=tuple(range(s,s+r))
    full=Psi(kap,(t,)*r,M)
    o,c=lead0(full,u)
    pred=(-1)**(s+r-1)*sp.Mul(*[w(k) for k in range(1,s)])*sp.Mul(*[w(m)**-2 for m in Mm])*Psi(kap[:-1],(t+1,)*(r-1),Mm)
    print(s,r,t,M,'order',o,'pred',t-s-r,'tail',-r,'ratio',sp.simplify(c/pred))
