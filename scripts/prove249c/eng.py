# exact engine: C(lam,mu,kap;c) via (D), Fractions
import itertools
from fractions import Fraction as Fr
def T0(kap,h,M,w):
    r=len(kap); x=[w[m] for m in h]; tails=[m for m in M if m not in h]
    e=Fr(1)
    for j in range(r):
        e*= (1-x[j]**-2)
        for k in range(1,kap[j]): e*=(1-w[k]/x[j])
    for j in range(r):
        for k in range(j+1,r): e*=(1-1/(x[j]*x[k]))/(1-x[j]/x[k])
    for t in tails:
        for j in range(r): e*=(1-1/(w[t]*x[j]))/(1-w[t]/x[j])
    return e
def Psi(kap,th,M,w):
    s=Fr(0)
    for h in itertools.permutations(M,len(kap)):
        if any(h[j]<kap[j] for j in range(len(kap))): continue
        x=[w[m] for m in h]; t=T0(kap,h,M,w)
        for j in range(len(kap)): t*=x[j]**th[j]
        s+=t
    return s
def F(kap,M,w):
    r=len(kap); p=Fr(1)
    for m in M: p*=w[m]**(2*r)
    return p*Psi(kap,[kap[j]-2*(j+1) for j in range(r)],M,w)
def cc(I,w):
    p=Fr(1)
    for i in I:
        p*=(-1)**i
        for k in range(1,i+1): p*=w[k]
    return p
def leq(A,B):  # top-aligned dominance A<=B
    A=sorted(A,reverse=True); B=sorted(B,reverse=True)
    if len(A)>len(B): return False
    return all(A[i]<=B[i] for i in range(len(A)))
def ell(kap,M,n,w):
    U=list(range(1,n+1)); kc=[i for i in U if i not in kap]; Mc=[i for i in U if i not in M]
    return cc(M,w)/cc(kap,w)*Psi(kc,[0]*len(kc),Mc,w)/(F(M,M,w)*Psi(Mc,[0]*len(Mc),Mc,w))
def subsets(n):
    for r in range(n+1):
        for S in itertools.combinations(range(1,n+1),r): yield list(S)
def Cfun(lam,mu,kap,c,n,w):
    s=Fr(0)
    for M in subsets(n):
        if not leq(M,kap) or not leq(lam,M) or not leq(mu,M): continue
        t=ell(kap,M,n,w)
        for m in M: t*=w[m]**c
        s+=t*F(lam,M,w)*F(mu,M,w)
    return s
