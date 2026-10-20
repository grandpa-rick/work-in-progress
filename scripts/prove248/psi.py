import sympy as sp, itertools
W=sp.symbols('w1:12')
def w(k): return W[k-1]
def T0(kap,h,M):
    r=len(kap); x=[w(m) for m in h]; tails=[m for m in M if m not in h]
    e=sp.Integer(1)
    for j in range(r):
        e*= (1-x[j]**-2)
        for k in range(1,kap[j]): e*=(1-w(k)/x[j])
    for j in range(r):
        for k in range(j+1,r): e*=(1-1/(x[j]*x[k]))/(1-x[j]/x[k])
    for t in tails:
        for j in range(r): e*=(1-1/(w(t)*x[j]))/(1-w(t)/x[j])
    return e
def Psi(kap,th,M):
    s=0
    for h in itertools.permutations(M,len(kap)):
        if any(h[j]<kap[j] for j in range(len(kap))): continue
        x=[w(m) for m in h]
        s+=T0(kap,h,M)*sp.Mul(*[x[j]**th[j] for j in range(len(kap))])
    return sp.factor(sp.cancel(sp.together(s)))
