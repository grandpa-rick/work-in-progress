# check (★): sigma Psi_lam(M+1) = z^2 Psi_lam(M) + z(1-z^2) w_M^{-2} Psi_{lam-}(M), lam ∋ 1, M ∌ 1, sigma: w1 -> 1/w1
import sys, random, itertools
from fractions import Fraction as Fr
sys.path.insert(0,'../prove249c')
from eng import Psi
random.seed(9); n=5
w={i:Fr(random.randint(2,40),random.randint(2,40)) for i in range(1,n+1)}
ws=dict(w); ws[1]=1/w[1]; z=w[1]
def tw(l): return [l[j]-2*(j+1) for j in range(len(l))]
ok=tot=0
for r in range(1,4):
  for lam in itertools.combinations(range(1,n+1),r):
    if 1 not in lam: continue
    lm=list(lam[1:])
    for k in range(0,n):
      for M in itertools.combinations(range(2,n+1),k):
        lhs=Psi(list(lam),tw(list(lam)),list(M)+[1],ws)
        wm=Fr(1)
        for m in M: wm*=w[m]
        rhs=z**2*Psi(list(lam),tw(list(lam)),list(M),w)+z*(1-z**2)*Psi(lm,tw(lm),list(M),w)/wm**2
        tot+=1; ok+=(lhs==rhs)
print('(★) exact:',ok,'/',tot)
# residue identity
T=[2,3,4]; zz=Fr(7,3)
S=Fr(0)
for i in T:
    t=(1-w[i]**-2)/(w[i]-1/zz)
    for u in T:
        if u!=i: t*=(w[i]-1/w[u])/(w[i]-w[u])
    S+=t
p1=Fr(1);p2=Fr(1)
for u in T: p1*=w[u]**-2; p2*=(w[u]-zz)/(w[u]*(1-zz*w[u]))
print('residue identity:', S== zz*p1-zz*p2)
