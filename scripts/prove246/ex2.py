import sys; sys.path.insert(0,'.')
from fast import *
N=4;a,b=1,2;nu=(1,2,4);lam=();mu=(3,)
for wt in (0,1):
  for sd in (3,11):
    B=Basis(N,a,b,sd,wt)
    terms=[]
    for M in [(1,2,3),(1,2,4)]:
        rho=B.coords({M:S.const(1)},nu).coef(-1)
        Mpp=tuple(c for c in M if c not in (a,b))
        Pv=(B.A(Mpp,lam)*B.A(Mpp,mu)).coef(0)
        terms.append((M,rho,Pv,rho*Pv))
    print(wt,sd,[(t[0],float(t[1]),float(t[2]),float(t[3])) for t in terms], float(terms[0][3]+terms[1][3]))
