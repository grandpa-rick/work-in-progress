from base import *
from check_conjT import antichain, comps
m=int(sys.argv[1]); KS,n,Ds=setup(m)
Ao=[[Fr(int(Ds[j]<=Ds[i] and antichain(Ds[i]-Ds[j]) and all(a!=b for a,b in Ds[i]-Ds[j]))) for j in range(n)] for i in range(n)]
Ae=[[Fr(int(Ds[j]<=Ds[i] and antichain(Ds[i]-Ds[j]))) for j in range(n)] for i in range(n)]
T=msqrt(mul(Ao,Ae,n),n)
sg=lambda i,j:(-1)**((sum(KS[i])-sum(KS[j]))%2)
Tt=[[sg(i,j)*T[i][j] for j in range(n)] for i in range(n)]
zo=[[Fr(int(Ds[j]<=Ds[i] and len(KS[i])==len(KS[j]))) for j in range(n)] for i in range(n)]
Q=mul(zo,Tt,n)
rim=lambda th: [x for x in th if (x[0]+1,x[1]+1) not in th]
def q1(i,j):
    if not Ds[j]<=Ds[i] or len(KS[i])!=len(KS[j])+1: return Fr(0)
    th=Ds[i]-Ds[j]
    if len(comps(th))!=1: return Fr(0)
    return -Fr(1,2**len(rim(th)))
bad=[(KS[i],KS[j],Q[i][j],q1(i,j)) for i in range(n) for j in range(n) if Q[i][j]-(i==j)!=q1(i,j)]
print(m,n,'Q formula failures',len(bad),bad[:5])
