"""Ask 1 (Clio 10-17): for an N=5 (c)-failure triple, FULL-q c^nu_{lam,mu} on every H_ab (all 10 pairs):
pole order in w along H_ab (numeric generic w, w_b = 1/(w_a(1-eps))), leading Laurent coefficient R(q),
its q-valuation at 0 and degree at q=inf.  Engine = day243/q488c/fullc.py (Motegi summand, E^C normalised)."""
import sys, os, itertools, sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../day243/q488c'))
from fullc import cfull, eps, q
idx=int(sys.argv[1]); N=5
T=eval(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'../../day243/q488c/n5fails_partial.txt')).read())
lam,mu,nu=T[idx]
def qval(r):
    n_,d_=sp.fraction(sp.cancel(r)); Pn=sp.Poly(n_,q); Pd=sp.Poly(d_,q)
    return min(x[0] for x in Pn.monoms())-min(x[0] for x in Pd.monoms()), Pn.degree()-Pd.degree()
for a,b in itertools.combinations(range(1,N+1),2):
    cv=sp.cancel(cfull(lam,mu,nu,N,a,b))
    n_,d_=sp.fraction(cv)
    on=min(m[0] for m in sp.Poly(n_,eps).monoms()); od=min(m[0] for m in sp.Poly(d_,eps).monoms())
    o=od-on
    if o<=0: print(lam,mu,nu,(a,b),'ord 0 (no pole)',flush=True); continue
    R=sp.cancel(cv*eps**o).subs(eps,0); v0,dinf=qval(R)
    extra=''
    if o>1:
        R2=sp.cancel((cv*eps**o-R)/eps).subs(eps,0)
        extra=' next-coeff qval=%s'%(qval(R2)[0] if R2!=0 else 'zero')
    print(lam,mu,nu,(a,b),'POLE ord',o,'lead qval0',v0,'deg_inf',dinf,extra,flush=True)
