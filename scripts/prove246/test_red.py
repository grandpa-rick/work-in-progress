import sys, random, itertools; sys.path.insert(0,'.')

import ptypes as ty
from fractions import Fraction as Fr
from fast import subsets, geq
R=random.Random(9)
def rv(n): return [Fr(R.randint(2,60),R.randint(2,60)) for _ in range(n)]
c=-6
# Test A4: n=b.  R^A_[n](lam,mu,nu;a,n)/B_[n-1](lam,mu,nu-;a) constant over (lam,mu)
n=4; a=2; vals=rv(n); vals[n-1]=1/vals[a-1]   # point on H (eps handled inside A)
nu=(1,2,3,4)
rat=set()
for lam,mu in [((),()),((),(1,)),((1,),(1,)),((),(2,)),((1,),(2,)),((),(3,)),((2,),(3,))]:
    r1=ty.A(n,vals,a,n,lam,mu,nu,c); r2=ty.B(n-1,vals[:n-1],a,lam,mu,nu[:-1],c)
    rat.add(r1/r2 if r2!=0 else ('B0',r1))
print('A4 ratios', [float(x) if not isinstance(x,tuple) else x for x in rat])
# Test B4: n=a. B_[n](lam,mu,kap;a=n)/C_[n-1](lam,mu,kap-; c-2)
n=4; vals=rv(n); kap=(1,3,4)
rat=set()
for lam,mu in [((),()),((),(1,)),((1,),(1,)),((),(2,)),((1,),(3,)),((),(3,))]:
    r1=ty.B(n,vals,n,lam,mu,kap,c); r2=ty.C(n-1,vals[:n-1],lam,mu,kap[:-1],c-2)
    rat.add(r1/r2 if r2!=0 else ('C0',r1))
print('B4 ratios', [float(x) if not isinstance(x,tuple) else x for x in rat])
