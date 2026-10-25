from base import *
m=int(sys.argv[1]); KS,n,Ds=setup(m)
L=lambda k: len(k)
Z=[[Fr(sum(1 for r in range(n) if Ds[j]<=Ds[r]<=Ds[i] and L(KS[r])==L(KS[j]))) if Ds[j]<=Ds[i] else Fr(0) for j in range(n)] for i in range(n)]
Y=msqrt(Z,n)
from collections import Counter
for i in range(n):
  for j in range(n):
    if Ds[j]<=Ds[i] and i!=j: print(KS[i],KS[j],Y[i][j], 'dl',L(KS[i])-L(KS[j]))
