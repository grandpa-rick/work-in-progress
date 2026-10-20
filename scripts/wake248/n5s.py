"""N=5 spot sample at c=7 (> 2N-4 = 6): lam=mu=empty (and mu=(1)), pair (1,2), several nu containing 5.
Same P1* test as cen.py (Weyl-stripped numerator one-signed)."""
import sys,time
from eng import *
N=5; C=setup(N)
from cen_core import analyse
s=Sys(N,1,2)
for (l,m,nu) in [((),(),(1,2,5)),((),(),(1,3,5)),((),(),(1,4,5)),((),(),(2,4,5)),((),(),(1,3,4,5)),((),(1,),(1,3,5)),((),(1,),(2,4,5))]:
    for c in (7,):
        if not EA(N,1,2,l,m,nu): print(l,m,nu,'not allowed'); continue
        t=time.time(); r=A(s,1,2,l,m,nu,c); d=analyse(r,N)
        print((1,2,l,m,nu),'c',c,{k:d[k] for k in d if k!='P'},'%.0fs'%(time.time()-t),flush=True)
