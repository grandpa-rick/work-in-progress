# validate order bookkeeping (head j: th_j-(r-1-j) [0-indexed], tail -r) against exact u->oo degree, own one-head/multi-head Psi
import sympy as sp, itertools, random
V=sp.Symbol('v')
def Psi(kap,th,M,x):   # primitive permutation sum (own implementation)
    s=0
    for h in itertools.permutations(M,len(kap)):
        if any(h[j]<kap[j] for j in range(len(kap))): continue
        X=[x[m] for m in h]; e=1
        for j in range(len(kap)):
            e*=(1-X[j]**-2)*X[j]**th[j]
            for k in range(1,kap[j]): e*=(1-x[k]/X[j])
            for k in range(j+1,len(kap)): e*=(1-1/(X[j]*X[k]))/(1-X[j]/X[k])
        for t in M:
            if t in h: continue
            for j in range(len(kap)): e*=(1-1/(x[t]*X[j]))/(1-x[t]/X[j])
        s+=e
    return s
def deg_inf(kap,th,M):
    random.seed(1); x={k:sp.Rational(random.randint(2,60),random.randint(2,60))+k for k in range(1,M[-1]+1)}
    x[M[-1]]=1/V
    num,den=sp.fraction(sp.cancel(sp.together(Psi(kap,th,M,x))))
    lo=lambda p: min(m[0] for m in sp.Poly(p,V).monoms())
    return lo(den)-lo(num)
for kap,th,M in [((2,),(-1,),(2,3,4)),((3,),(-1,),(3,4,5)),((1,),(-1,),(1,2,3)),((2,3),(-2,-2),(2,3,4,5)),
                 ((1,3),(-1,-2),(1,2,3,4)),((1,2),(-2,-2),(1,2,3,4)),((3,4),(-1,-1),(3,4,5,6)),((2,3),(-1,-1),(2,3,4,5)),((2,4),(-2,-3),(2,3,4,5))]:
    r=len(kap); heads=[th[j]-(r-1-j) for j in range(r)]; pred=max(max(heads),-r)
    tie=max(heads)==-r
    print(kap,th,M,'exact deg',deg_inf(kap,th,M),'pred',pred,'TIE' if tie else '')
