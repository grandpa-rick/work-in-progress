import sympy as sp, random
from q532 import Psi, pt
def R(x): return {k:sp.Rational(v.numerator,v.denominator) for k,v in x.items()}
def chk(r,S_):
  ok=[]
  for sd in (3,4):
    random.seed(sd); x=R(pt(max(S_)+1)); e=sp.Mul(*[x[m] for m in S_])
    L=(Psi(tuple(range(1,r)),(-r,)*(r-1),S_,x) if r>1 else 1)+(-1)**r*Psi(tuple(range(1,r+1)),(1-r,)*r,S_,x)
    Rr=e**-2*(Psi(tuple(range(1,r)),(2-r,)*(r-1),S_,x) if r>1 else 1)
    ok.append(L==Rr)
  return ok
for r,S_ in [(1,(1,2,3)),(2,(1,2,3,4)),(3,(1,2,3,4)),(3,(1,2,4,5)),(4,(1,2,3,4,5))]:
  print('P1-infinity conj r',r,'S',S_,chk(r,S_))
# breaking control: wrong twist on RHS
random.seed(5); x=R(pt(6)); S_=(1,2,3,4); e=sp.Mul(*[x[m] for m in S_])
L=Psi((1,2),(-3,-3),S_,x)-Psi((1,2,3),(-2,-2,-2),S_,x)
print('control RHS twist 3-r (should be False):', L==e**-2*Psi((1,2),(0,0),S_,x))
# s=2 recurring tie (2),(-1): head 1, tail -Psi_{(2),0}; fit?
for S_ in [(2,3,4),(1,2,3,4)]:
  random.seed(6); x=R(pt(5)); e=sp.Mul(*[x[m] for m in S_])
  T=1-Psi((2,),(0,),S_,x)
  fits=[(kap,t,k) for kap in [(1,),(2,)] for t in range(-3,4) for k in range(-4,3) if Psi(kap,(t,),S_,x)!=0 and T==e**k*Psi(kap,(t,),S_,x)]
  print('s=2 tie 1-Psi_(2),0 on',S_,'single fits:',fits, ' two-term form e^-2 + w1*Psi_(1),-1 holds:', T==e**-2+x[1]*Psi((1,),(-1,),S_,x))
