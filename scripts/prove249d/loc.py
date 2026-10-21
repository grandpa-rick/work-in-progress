# Termwise Theorem D: rho_{kappa L} == A_kappa * a^kappa_L(y) * w_L^2, w_m = y_{n+1-m}
# a^kappa_L = localization coefficients of D_{w_kappa}(e^{t omega}), computed with coefficients stored on W-orbit of y
import sys, random, itertools
from fractions import Fraction as Fr
sys.path.insert(0,'../prove249c')
from eng import ell, leq, subsets
from dem import word
n=int(sys.argv[1]); random.seed(int(sys.argv[2]) if len(sys.argv)>2 else 7)
y0=tuple(Fr(random.randint(2,50),random.randint(2,50)) for _ in range(n))
# Weyl group C_n acting on points y (multiplicative): signed permutations: y -> (y_{p(i)}^{s_i})
def act(i,y):
    y=list(y)
    if i<n: y[i-1],y[i]=y[i],y[i-1]
    else: y[n-1]=1/y[n-1]
    return tuple(y)
orbit={y0}; fr=[y0]
while fr:
    nf=[]
    for p in fr:
        for i in range(1,n+1):
            q=act(i,p)
            if q not in orbit: orbit.add(q); nf.append(q)
    fr=nf
orbit=list(orbit)
def ema(i,y):  # e^{-alpha_i} at point y
    return y[i]/y[i-1] if i<n else y[n-1]**-2
def sv(i,v):
    v=list(v)
    if i<n: v[i-1],v[i]=v[i],v[i-1]
    else: v[n-1]=-v[n-1]
    return tuple(v)
def Dop(cls,i):
    # cls: v -> {point: value}; result coefficient at v': (a_{v'}(y) - e^{-a}(y) a_{s v'}(s y))/(1-e^{-a}(y))
    keys=set(cls)|{sv(i,v) for v in cls}
    out={}
    for vp in keys:
        a=cls.get(vp); b=cls.get(sv(i,vp))
        d={}
        for y in orbit:
            e=ema(i,y)
            val=(a[y] if a else 0) - e*(b[act(i,y)] if b else 0)
            d[y]=val/(1-e)
        if any(d.values()): out[vp]=d
    return out
bad=0; tot=0
for kap in subsets(n):
    cls={tuple([1]*n):{y:Fr(1) for y in orbit}}
    for i in reversed(word(set(kap),n)): cls=Dop(cls,i)
    w={m:y0[n-m] for m in range(1,n+1)}
    A=Fr((-1)**sum(kap))
    for k in kap:
        for j in range(1,k+1): A/=w[j]
    for L in subsets(n):
        v=tuple(-1 if (n+1-p) in L else 1 for p in range(1,n+1))
        aL=cls[v][y0] if v in cls else Fr(0)
        wL2=Fr(1)
        for m in L: wL2*=w[m]**2
        r=ell(kap,L,n,w) if leq(L,kap) else Fr(0)
        tot+=1
        if r!=A*aL*wL2: bad+=1; print('MISMATCH',kap,L,r,A*aL*wL2)
print('n',n,'pairs (kappa,L)',tot,'mismatches',bad)
