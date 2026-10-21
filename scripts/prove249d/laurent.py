import sympy as sp, itertools, sys
n=int(sys.argv[1]); W=sp.symbols('w1:%d'%(n+1)); w={i+1:W[i] for i in range(n)}
def T0(kap,h,M):
    r=len(kap); x=[w[m] for m in h]; tails=[m for m in M if m not in h]; e=1
    for j in range(r):
        e*=(1-x[j]**-2)
        for k in range(1,kap[j]): e*=(1-w[k]/x[j])
    for j in range(r):
        for k in range(j+1,r): e*=(1-1/(x[j]*x[k]))/(1-x[j]/x[k])
    for t in tails:
        for j in range(r): e*=(1-1/(w[t]*x[j]))/(1-w[t]/x[j])
    return e
def F(kap,M):
    r=len(kap); s=0
    for h in itertools.permutations(M,r):
        if any(h[j]<kap[j] for j in range(r)): continue
        t=T0(kap,h,M)
        for j in range(r): t*=w[h[j]]**(kap[j]-2*(j+1))
        s+=t
    p=1
    for m in M: p*=w[m]**(2*r)
    return p*s
bad=0;tot=0
for r in range(1,n+1):
  for kap in itertools.combinations(range(1,n+1),r):
    for rm in range(r,n+1):
      for M in itertools.combinations(range(1,n+1),rm):
        v=sp.factor(sp.together(F(list(kap),list(M))))
        if v==0: continue
        num,den=sp.fraction(v); tot+=1
        if len(sp.Poly(den,*W).terms())!=1: bad+=1; print('NON-LAURENT',kap,M,den)
print('n',n,'nonzero F values',tot,'non-Laurent',bad)
