import sys, random; sys.path.insert(0,'.')
import ptypes as ty
from fractions import Fraction as Fr
R=random.Random(6)
def rv(n): return [Fr(R.randint(2,60),R.randint(2,60)) for _ in range(n)]
c=-6
def Kl(l,vals):
    if not l: return Fr(1)
    r=Fr((-1)**l[-1])
    for k in range(1,l[-1]): r*=vals[k-1]
    return r
mn=lambda t: t[:-1]
n=4; base=rv(n)
print('B2 (a=2, w4->0): B*u^-(c-delta)/(K K B_[3])')
for lam,mu,kap in [((),(),(1,2,3,4)),((),(3,),(2,3,4)),((1,),(1,),(1,2,3,4)),((),(4,),(1,2,4))]:
    out=[]
    for k in (8,14):
        v=list(base); u=Fr(1,10**k); v[3]=u; d=(1 if lam else 0)+(1 if mu else 0)
        r=ty.B(n,v,2,lam,mu,kap,c)*u**(-(c-d)); rp=ty.B(n-1,v[:3],2,mn(lam),mn(mu),mn(kap),c)
        out.append(float(r/(Kl(lam,v)*Kl(mu,v)*rp)))
    print('  ',lam,mu,kap,out)
print('B3 (a=2, stuck, w4->inf): B*u/B_[3](kap~)')
for lam,mu,kap,kt in [((),(),(3,4),(2,3)),((),(1,),(1,3,4),(1,2,3))]:
    out=[]
    for k in (8,14):
        v=list(base); u=Fr(10**k); v[3]=u
        r=ty.B(n,v,2,lam,mu,kap,c)*u; rp=ty.B(n-1,v[:3],2,lam,mu,kt,c)
        out.append(float(r/rp) if rp!=0 else ('rp0',float(r)))
    print('  ',lam,mu,kap,'->',kt,out)
print('C2 (w4->0): C*u^-(c+1-delta)/(K K C_[3])')
for lam,mu,kap in [((),(),(1,4)),((1,),(2,),(1,2,4)),((),(4,),(3,4))]:
    out=[]
    for k in (8,14):
        v=list(base); u=Fr(1,10**k); v[3]=u; d=(1 if lam else 0)+(1 if mu else 0)
        r=ty.C(n,v,lam,mu,kap,c)*u**(-(c+1-d)); rp=ty.C(n-1,v[:3],mn(lam),mn(mu),mn(kap),c)
        out.append(float(r/(Kl(lam,v)*Kl(mu,v)*rp)))
    print('  ',lam,mu,kap,out)
print('C3 (stuck: lam=mu=(), E={()}), w4->inf: C*u/C_[3](kap~)')
for kap,kt in [((4,),(3,)),((2,4),(2,3)),((3,4),(2,3))]:
    out=[]
    for k in (8,14):
        v=list(base); u=Fr(10**k); v[3]=u
        out.append(float(ty.C(n,v,(),(),kap,c)*u/ty.C(n-1,v[:3],(),(),kt,c)))
    print('  ',kap,'->',kt,out)
