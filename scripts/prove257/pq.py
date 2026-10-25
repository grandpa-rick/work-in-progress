from base import *
m=int(sys.argv[1]); KS,n,Ds=setup(m)
M=[[Fr(0)]*n for _ in range(n)]
Ao=[[Fr(0)]*n for _ in range(n)]
from check_conjT import antichain
Ao=[[Fr(int(Ds[j]<=Ds[i] and antichain(Ds[i]-Ds[j]) and all(a!=b for a,b in Ds[i]-Ds[j]))) for j in range(n)] for i in range(n)]
Ae=[[Fr(int(Ds[j]<=Ds[i] and antichain(Ds[i]-Ds[j]))) for j in range(n)] for i in range(n)]
T=msqrt(mul(Ao,Ae,n),n)
sg=lambda i,j:(-1)**(sum(KS[i])-sum(KS[j]))
Tt=[[sg(i,j)*T[i][j] for j in range(n)] for i in range(n)]
z=[[Fr(int(Ds[j]<=Ds[i])) for j in range(n)] for i in range(n)]
zo=[[Fr(int(Ds[j]<=Ds[i] and len(KS[i])==len(KS[j]))) for j in range(n)] for i in range(n)]
P=mul(Tt,z,n); Q=mul(zo,Tt,n)
which=sys.argv[2]
X=P if which=='P' else Q
for i in range(n):
  for j in range(n):
    if X[i][j]!=0 and i!=j: print(KS[i],KS[j],X[i][j])
