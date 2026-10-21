# reflection probe: Psi_{[r],t}(S) =? +-e_S^k Psi_{[r],c-t}(S)
import sympy as sp, random
from q532 import Psi, pt
def R(x): return {k:sp.Rational(v.numerator,v.denominator) for k,v in x.items()}
for r,S in [(1,(1,2,3)),(2,(1,2,3,4)),(3,(1,2,3,4,5))]:
  res=[]
  for t in range(-2,3):
    found=None
    for c in [c for c in range(-2*r-4,5) if c!=2*t]:
      ok=True; kk=None
      for sd in (11,12):
        random.seed(sd); x=R(pt(max(S))); e=sp.Mul(*[x[m] for m in S])
        A=Psi(tuple(range(1,r+1)),(t,)*r,S,x); B=Psi(tuple(range(1,r+1)),(c-t,)*r,S,x)
        if B==0: ok=False;break
        q=A/B; m=[k for k in range(-12,13) if q in (e**k,-e**k)]
        if not m or (kk is not None and m[0]!=kk): ok=False;break
        kk=m[0]
      if ok: found=(c,kk);break
    res.append((t,found))
  print('r',r,res)
