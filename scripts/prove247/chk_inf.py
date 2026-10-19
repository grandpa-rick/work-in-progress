# new step: mixed case, lam=mu=empty, c=-1: u*X_n(nu) -> const(m*) * X_{n-1}(nu~)  as u->inf
import sys; sys.path.insert(0,'../prove246')
from fractions import Fraction as Fr
import random, ptypes as ty
R=random.Random(4)
c=-1
def tilde(nu,n):
    comp=[m for m in range(1,n+1) if m not in nu]
    if not comp: return tuple(m for m in nu if m!=n), None
    ms=max(comp); return tuple(sorted([m for m in nu if m!=n]+[ms])), ms
for n,a,b in [(4,1,2),(5,1,2),(5,2,3)]:
  base=[Fr(R.randint(2,40),R.randint(2,40)) for _ in range(n-1)]
  for U in (Fr(10**6),Fr(10**9)):
    vals=base+[U]; res={}
    for nu in [(1,2,4),(1,2,3,n),(1,2,n),(2,3,n),(1,2,3,4,5)[:n], (1,3,n),(2,n)]:
        if n not in nu or a not in nu and b not in nu: pass
        nt,ms=tilde(nu,n)
        if not (ty.geq(nu,(a,b)) if hasattr(ty,'geq') else True): continue
        try:
            r1=ty.A(n,vals,a,b,(),(),nu,c); r2=ty.A(n-1,base,a,b,(),(),nt,c)
        except Exception as e: continue
        if r2==0: res[nu]=('r2=0',float(r1)); continue
        res.setdefault(ms,[]).append((nu,round(float(U*r1/r2),5)))
    print(n,a,b,float(U),res)
