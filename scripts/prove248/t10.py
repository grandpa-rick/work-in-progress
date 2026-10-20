# Lemma P* check, numeric (random rational point) to be fast
from psi import *
import random
random.seed(5)
vals={w(k):sp.Rational(random.randint(2,40),random.randint(2,40)) for k in range(1,9)}
def Psin(kap,th,M,uval=None):
    s=0
    for h in itertools.permutations(M,len(kap)):
        if any(h[j]<kap[j] for j in range(len(kap))): continue
        s+=(T0(kap,h,M)*sp.Mul(*[w(h[j])**th[j] for j in range(len(kap))]))
    return s
def lead_num(kap,th,M,o):
    n=M[-1]; u=w(n); v={k:v for k,v in vals.items() if k!=u}
    e=sp.cancel(sp.together(Psin(kap,th,M).subs(v)))
    f=sp.cancel(e*u**(-o+1)); f0=f.subs(u,0)
    f1=sp.cancel(f/u).subs(u,0) if f0==0 else None
    return f1, f0
for kap,th,M in [((2,3,6),(0,0,4),(2,3,4,6,7)),((1,2,4),(0,0,2),(1,2,3,4,5))]:
    n=M[-1]; Mm=M[:-1]; r=len(kap)
    tb=th[0]; o=tb-kap[0]-1-(r-1)
    nl=(kap[0],kap[2]); nt=(tb+1,th[2]+1)
    pred=((-1)**(kap[0]+r-1)*sp.Mul(*[w(k) for k in range(1,kap[0])])*sp.Mul(*[w(m)**-2 for m in Mm])*Psin(nl,nt,Mm)).subs(vals)
    c,cbelow=lead_num(kap,th,M,o)
    print(kap,th,M,'order',o,'coef-pred',sp.simplify(c-pred),'u^(o-1) coef (should be 0)',cbelow, 'other orders: top head',th[2]-kap[2]-1,'tail',-r)
