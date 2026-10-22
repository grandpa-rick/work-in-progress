import sympy as sp
from robust import calH, K
tau=sp.Symbol('tau')
ok=tot=0
for lam in [[3,1],[4,2],[4,1],[5,2]]:
    k1=lam[0]; lp=lam[1:]
    r=sp.Rational((-1)**k1,2**(k1+1))
    for c in [sp.Rational(n,2) for n in range(2,k1+2)]:
        if 2*c>k1+1: continue
        for T in ([],[tau]):
            lhs=calH(lam,c,T)
            rhs=-r*(1-(-1)**(2*c+len(T)))*calH(lp,c-sp.Rational(1,2),T)
            for tt in T:
                rest=[u for u in T if u!=tt]
                rhs+= -(-1)**k1*tt**(k1+1-2*c)*(tt-1)**(-k1)*sp.Mul(*[K(u,1/tt) for u in rest])*calH(lp,c-sp.Rational(1,2),rest)
            tot+=1; good=sp.simplify(lhs-rhs)==0; ok+=good
            if not good: print('MISS',lam,c,T)
print('Prop M',ok,'/',tot)
