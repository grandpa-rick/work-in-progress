from half import *
from fractions import Fraction as Fr
import sys
N=int(sys.argv[1])
idx={k:i for i,k in enumerate(KS)}
n=len(KS)
I=[[Fr(int(i==j)) for j in range(n)] for i in range(n)]
C=[[Fr(0)]*n for _ in range(n)]
for k,d in SUP.items():
    for r,v in d.items(): C[idx[k]][idx[r]]=Fr(v)
def mul(A,B): return [[sum(A[i][l]*B[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
# T = sum binom(1/2,k) C^k
T=[row[:] for row in I]; P=[row[:] for row in I]; b=Fr(1)
for k in range(1,N*N+2):
    P=mul(P,C); b=b*(Fr(1,2)-(k-1))/k
    T=[[T[i][j]+b*P[i][j] for j in range(n)] for i in range(n)]
for k in KS:
    if k: print(k, {(r if r else '0'):str(T[idx[k]][idx[r]]) for r in KS if T[idx[k]][idx[r]]!=0})
# check anchor
for k in KS:
    v=sum(T[idx[k]][idx[r]]*(-1)**sum(r) for r in KS)
    assert v==Fr((-1)**sum(k),2**(k[0] if k else 0)), (k,v)
print('anchor identity T s_0 = s_1/2 OK')
