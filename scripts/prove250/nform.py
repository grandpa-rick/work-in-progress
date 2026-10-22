import sys, itertools, sympy as sp
n=int(sys.argv[1])
W=sp.symbols('w1:%d'%(n+1)); w={m:W[m-1] for m in range(1,n+1)}
def T(kap,th,M):  # single term h=identity on M=kap
    e=1; r=len(kap); x=[w[m] for m in kap]
    for j in range(r):
        e*=(1-x[j]**-2)*x[j]**th[j]
        for k in range(1,kap[j]): e*=(1-w[k]/x[j])
    for j in range(r):
        for k in range(j+1,r): e*=(1-1/(x[j]*x[k]))/(1-x[j]/x[k])
    return e
def cc(I):
    p=1
    for i in I:
        p*=(-1)**i
        for k in range(1,i+1): p*=w[k]
    return p
for r in range(n+1):
  for L in itertools.combinations(range(1,n+1),r):
    L=list(L); C=[i for i in range(1,n+1) if i not in L]
    u={k:(1/w[k] if k in L else w[k]) for k in w}
    N=cc(L)/( sp.prod([w[m]**(2*r+2) for m in L]) * T(L,[L[j]-2*(j+1) for j in range(r)],L) * T(C,[0]*len(C),C))
    E=sp.prod([1-u[k]*u[l] for k in range(1,n+1) for l in range(k,n+1)])
    E2=sp.prod([1-u[k]/u[l] for k in range(1,n+1) for l in range(k+1,n+1)])
    print(L, sp.factor(sp.simplify(N*E)), '|', sp.factor(sp.simplify(N*E*E2)))
