import sys, random; sys.path.insert(0,'.')
import ptypes as ty
from fractions import Fraction as Fr
R=random.Random(4)
def rv(n): return [Fr(R.randint(2,60),R.randint(2,60)) for _ in range(n)]
c=-6
def Kl(l,vals):
    if not l: return Fr(1)
    r=Fr((-1)**l[-1])
    for k in range(1,l[-1]): r*=vals[k-1]
    return r
mn=lambda t: t[:-1]
# A2: n=4 top, pair (1,2), nu containing 4, in-terms exist
n=4; base=rv(n); base[1]=1/base[0]
print('A2 (w4->0): R*u^-(c+1-delta) / (K_lam K_mu R_[3](lam-,mu-,nu-))')
for lam,mu,nu in [((),(),(1,2,3,4)),((),(3,),(1,2,4)),((1,),(3,),(1,2,3,4)),((),(4,),(1,2,4)),((2,),(2,),(1,2,3,4))]:
    out=[]
    for k in (8,14):
        v=list(base); u=Fr(1,10**k); v[3]=u
        d=(1 if lam else 0)+(1 if mu else 0)
        r=ty.A(n,v,1,2,lam,mu,nu,c)*u**(-(c+1-d))
        rp=ty.A(n-1,v[:3],1,2,mn(lam),mn(mu),mn(nu),c)
        out.append(float(r/(Kl(lam,v)*Kl(mu,v)*rp)))
    print('  ',lam,mu,nu,out)
# A3: stuck case nu=(3,4), pair (1,2): E={()}; nu~ = (3) u {m*}: nu^c=(1,2) -> m*=2 -> nu~=(2,3)
print('A3 (w4->inf): R*u / R_[3](nu~)')
for lam,mu,nu,nut in [((),(),(3,4),(2,3)),((),(),(1,3,4),(1,2,3)),((),(1,),(2,3,4),(1,2,3))]:
    out=[]
    for k in (8,14):
        v=list(base); u=Fr(10**k); v[3]=u
        r=ty.A(n,v,1,2,lam,mu,nu,c)*u
        rp=ty.A(n-1,v[:3],1,2,lam,mu,nut,c)
        out.append(float(r/rp) if rp!=0 else ('rp0',float(r)))
    print('  ',lam,mu,nu,'->',nut,out)
