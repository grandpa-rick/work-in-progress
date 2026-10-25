import itertools, sys
from fractions import Fraction as Fr
def strict(m):
    return [tuple(sorted(c, reverse=True)) for L in range(m+1) for c in itertools.combinations(range(1, m+1), L)]
def D(k): return frozenset((a, b) for a, ka in enumerate(k) for b in range(a, a+ka))
def setup(m):
    KS = sorted(strict(m), key=sum); n=len(KS); Ds=[D(k) for k in KS]
    return KS,n,Ds
def mul(A,B,n): return [[sum(A[i][l]*B[l][j] for l in range(n) if A[i][l]) for j in range(n)] for i in range(n)]
def msqrt(M,n,p=Fr(1,2)):
    C=[[M[i][j]-(i==j) for j in range(n)] for i in range(n)]
    T=[[Fr(int(i==j)) for j in range(n)] for i in range(n)]; P=[r[:] for r in T]; b=Fr(1)
    for k in range(1,n+1):
        P=mul(P,C,n); b=b*(p-(k-1))/k
        if all(x==0 for r in P for x in r): break
        T=[[T[i][j]+b*P[i][j] for j in range(n)] for i in range(n)]
    return T
