# lambda=[r] (s=1) tie at u->oo: is head+tail a single instance (lambda-side P1 analogue)? fit e_S^k * Psi_{kap',t}(S)
import sympy as sp, random
from q532 import Psi, pt
def R(x): return {k:sp.Rational(v.numerator,v.denominator) for k,v in x.items()}
for r,M in [(2,(1,2,3,4)),(2,(1,2,3,4,5)),(3,(1,2,3,4,5))]:
  hits=[]
  pts=[]
  for sd in (1,2):
    random.seed(sd); x=R(pt(max(M))); Mm=M[:-1]
    S=Psi(tuple(range(1,r)),(-r,)*(r-1),Mm,x)+(-1)**r*Psi(tuple(range(1,r+1)),(1-r,)*r,Mm,x)
    eS=sp.Mul(*[x[m] for m in Mm]); pts.append((x,S,eS))
  for kap in [tuple(range(1,r)),tuple(range(1,r+1)),tuple(range(2,r+1)),tuple(range(1,r-1))]:
    for t in range(-r-2,3):
      vals=[]
      for x,S,eS in pts:
        P=Psi(kap,(t,)*len(kap),M[:-1],x) if kap else sp.Integer(1)
        if P==0: vals=None;break
        q=S/P
        k=None
        for kk in range(-6,7):
          if q==eS**kk or q==-eS**kk: k=(kk,sp.sign(q/eS**kk)); break
        vals.append(k)
      if vals and vals[0] is not None and vals[0]==vals[1]: hits.append((kap,t,vals[0]))
  print('r',r,'M',M,'single-instance fits (kap,t,(k,sign)):',hits)
