# Day 213: DS length-2 from the 207b Pieri formula. Exact sympy only.
from sympy import symbols, Rational, prod, cancel, simplify, factor, Integer
s,t,w=symbols('s t w')
def qp(a,n): return prod([1-a*t**i for i in range(n)]) if n>0 else Integer(1)
def alpha(j):
    if j<0: return Integer(0)
    return prod([s-t**i for i in range(1,j+1)])/qp(t,j) if j>0 else Integer(1)
def c(n,j):
    if j<0 or j>n: return Integer(0)
    return qp(s,n-j)/qp(t,n-j)*(alpha(j)-s*t**(n-j)*alpha(j-1))
def F(n,x): return sum(c(n,j)*x**j for j in range(n+1))
ok=True
# (V) F_n|_{s=1} = 0 identically in w, n>=1 ; F_0 = 1
for n in range(0,8):
    v=cancel(F(n,w).subs(s,1))
    if (n>=1 and v!=0) or (n==0 and v!=1): ok=False; print("FAIL V",n,v)
# per (k,r) with k<=r: coefficients of e_(r+k-b, b)
for k in range(1,5):
    for r in range(k,6):
        coeffs={}
        for b in range(k+1):
            mu=(r+k-b,b)
            assert mu not in coeffs
            coeffs[mu]=cancel(s**b*F(k-b,t**(r-b)))
        # lead
        if cancel(coeffs[(r,k)]-s**k)!=0: ok=False; print("FAIL lead",k,r)
        for mu,cf in coeffs.items():
            if mu==(r,k): continue
            if cancel(cf.subs(s,1))!=0: ok=False; print("FAIL q=1",k,r,mu)
            if cf==0 or cancel(cf.subs(t,0)-(1-s)*s**mu[1])!=0: ok=False; print("FAIL nonzero/t0",k,r,mu)
# swap symmetry of collected RHS, small cases (sanity, already proved in 207b s8)
from collections import defaultdict
def rhs(k,r):
    d=defaultdict(lambda: Integer(0))
    for b in range(k+1):
        mu=tuple(sorted((b,r+k-b),reverse=True)); mu=tuple(x for x in mu if x>0)
        d[mu]+=s**b*F(k-b,t**(r-b))
    return {m:cancel(v) for m,v in d.items() if cancel(v)!=0}
for k in range(1,5):
    for r in range(1,k):
        A=rhs(k,r); B=rhs(r,k)
        if set(A)!=set(B) or any(cancel(A[m]-B[m])!=0 for m in A): ok=False; print("FAIL swap",k,r)
print("ALL OK" if ok else "FAILURES")
