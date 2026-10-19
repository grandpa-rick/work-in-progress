import sys; sys.path.insert(0,'../prove246')
from fractions import Fraction as Fr
import random, ptypes as ty
R=random.Random(5)
for n,a,b,c in [(4,1,2,5),(4,1,3,6),(5,1,2,7),(5,2,4,9)]:
  base=[Fr(R.randint(2,40),R.randint(2,40)) for _ in range(n-1)]
  for U in (Fr(10**6),Fr(10**9)):
    out=[]
    for nu in [(1,2,3,4,5)[:n],(1,2,n),(2,3,n),(a,b,n),(1,3,4,5)[:n-1]+(n,) ]:
        nu=tuple(sorted(set(nu)))
        if not (len(nu)>=2): continue
        try: r1=ty.A(n,base+[U],a,b,(),(),nu,c); r2=ty.A(n-1,base,a,b,(),(),nu[:-1],c-2)
        except Exception as e: out.append((nu,'err')); continue
        out.append((nu, round(float(r1/U**(c-n)/r2),6) if r2 else ('r2=0',float(r1))))
    print(n,a,b,c,float(U),out)
