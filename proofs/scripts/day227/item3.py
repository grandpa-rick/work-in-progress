# Item 3: lin_e P_lam = (-1)^{n-1}(1-t^n) X^lam_{(n)}/b_lam  [green.py, characters + Kostka-Foulkes]
# vs one-string residue (Day 225 Cor 2.3): (-1)^{n-k}[n]/[k] P_rho(1,t,..,t^{k-1}); P_rho by brute symmetrization
# at the principal point (denominators z_i-z_j nonzero there), exact Fractions, t in {2/3, 3, -5/2}.
import sys; sys.path.insert(0,'../day226')
import sympy as sp, itertools
from fractions import Fraction as Fr
from collections import Counter
from green import X as Xg, parts, t, b
def Pprin(rho,k,T):
    z=[T**i for i in range(k)]; r=list(rho)+[0]*(k-len(rho))
    v=Fr(1)
    for m in Counter(r).values():
        for i in range(1,m+1): v*=sum(T**j for j in range(i))
    s=Fr(0)
    for w in itertools.permutations(range(k)):
        zz=[z[w[i]] for i in range(k)]; term=Fr(1)
        for i in range(k):
            term*=zz[i]**r[i]
            for j in range(i+1,k): term*=(zz[i]-T*zz[j])/(zz[i]-zz[j])
        s+=term
    return s/v
q=lambda m,T: sum(T**i for i in range(m))
NMAX=int(sys.argv[1]) if len(sys.argv)>1 else 7
ok=0;tot=0
for n in range(1,NMAX+1):
    for lam in parts(n):
        L=sp.cancel((-1)**(n-1)*(1-t**n)*Xg(lam,(n,))/b(lam))
        k=len(lam); rho=tuple(x-1 for x in lam if x>1)
        for T in [Fr(2,3),Fr(3),Fr(-5,2)]:
            v=L.subs(t,sp.Rational(T.numerator,T.denominator)); lhs=Fr(int(v.p),int(v.q))
            rhs=(-1)**(n-k)*q(n,T)/q(k,T)*Pprin(rho,k,T)
            tot+=1; ok+= lhs==rhs
print('item3',ok,'/',tot)
