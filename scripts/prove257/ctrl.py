# Verify the proof's objects: X (seq formula), P (seq formula), (I1) XP=I, (I2) X zeta = zeta_o P,
# T' = S zeta_o^{-1} X S equals T = sqrt(A_o A_e) (binomial series) and the printed rule; (a).
from base import *
from check_conjT import antichain, rule
m=int(sys.argv[1]); KS,n,Ds=setup(m)
part=lambda k,r: k[r-1] if r<=len(k) else 0
def sub(nu,ka): return all(part(nu,r)<=part(ka,r) for r in range(1,max(len(nu),len(ka))+1))
def r0(ka,nu):
    r=1
    while part(nu,r)==part(ka,r): r+=1
    return r
def Xf(ka,nu):
    if ka==nu: return Fr(1)
    if len(nu)!=len(ka)-1: return Fr(0)
    r=r0(ka,nu)
    if r>len(ka): return Fr(0)
    if not all(part(nu,s)<=part(ka,s+1) for s in range(r,len(ka))): return Fr(0)
    return -Fr(1,2**part(ka,r))*(Fr(3,2) if len(ka)>=2 else 1)
def Pf(ka,nu):
    if ka==nu: return Fr(1)
    if len(nu)>=len(ka): return Fr(0)
    r=r0(ka,nu)
    if not all(part(nu,s)<=part(ka,s+1) for s in range(r,len(ka)+1)): return Fr(0)
    return Fr(1,2**part(ka,r))
X=[[Xf(KS[i],KS[j]) for j in range(n)] for i in range(n)]
P=[[Pf(KS[i],KS[j]) for j in range(n)] for i in range(n)]
z=[[Fr(int(sub(KS[j],KS[i]))) for j in range(n)] for i in range(n)]
assert all(z[i][j]==int(Ds[j]<=Ds[i]) for i in range(n) for j in range(n))  # containment = partwise
zo=[[z[i][j]*(len(KS[i])==len(KS[j])) for j in range(n)] for i in range(n)]
I=[[Fr(int(i==j)) for j in range(n)] for i in range(n)]
print('I1 XP=I', mul(X,P,n)==I)
print('I2 X z = zo P', mul(X,z,n)==mul(zo,P,n))
Ao=[[Fr(int(Ds[j]<=Ds[i] and antichain(Ds[i]-Ds[j]) and all(a!=b for a,b in Ds[i]-Ds[j]))) for j in range(n)] for i in range(n)]
Ae=[[Fr(int(Ds[j]<=Ds[i] and antichain(Ds[i]-Ds[j]))) for j in range(n)] for i in range(n)]
sg=[(-1)**(sum(k)%2) for k in KS]
S=lambda M:[[sg[i]*sg[j]*M[i][j] for j in range(n)] for i in range(n)]
# Moebius identifications
print('Ao zo = I (signed)', mul(S(Ao),zo,n)==I, ' Ae z = I (signed)', mul(S(Ae),z,n)==I)
Tp=mul(Ao,S(X),n)   # S zo^{-1} X S = (S zo^{-1} S)(S X S) = Ao * S(X)
print('Tp^2 = AoAe', mul(Tp,Tp,n)==mul(Ao,Ae,n))
print('Tp = rule', all(Tp[i][j]==(rule(Ds[i]-Ds[j],'NEW') if Ds[j]<=Ds[i] else 0) for i in range(n) for j in range(n)))
if m<=6:
    T=msqrt(mul(Ao,Ae,n),n); print('Tp = binomial-series T', Tp==T)
s0=[Fr(s) for s in sg]; sh=[Fr(sg[i],2**(KS[i][0] if KS[i] else 0)) for i in range(n)]
print('(a) Tp s0 = s_1/2', [sum(Tp[i][j]*s0[j] for j in range(n)) for i in range(n)]==sh)
