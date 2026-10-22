import sympy as sp, itertools
from robust import K
def term(lam,c,S):
    l=len(lam); t=sp.symbols('t1:%d'%(l+1)); expr=1
    for s in range(1,l+1):
        ts=t[s-1]; m=ts/((ts**2-1)*(ts-1)**lam[s-1])
        if s in S:
            rho=ts**(2*c-1-s)*sp.Mul(*[K(t[r-1],ts) for r in range(1,s)])
            expr*=-m*rho
        else:
            # keep m but kernels into later S-variables are inside their rho already
            expr*=m
    for s in range(l,0,-1): expr=sp.residue(sp.together(expr),t[s-1],1)
    return sp.simplify(expr)
lam=[3,2,1]; c=sp.Rational(5,2)
full=0
for q in range(4):
    for S in itertools.combinations([1,2,3],q):
        v=term(lam,c,S); full+=v
        # decoupled prediction: rows not in S -> factor -r_k, and kernels to them removed
        print(S, v)
print('sum',sp.simplify(full))
def JS(lam,c,S):
    t={s:sp.Symbol('t%d'%s) for s in S}; expr=1
    for s in S:
        ts=t[s]; m=ts/((ts**2-1)*(ts-1)**lam[s-1])
        expr*=m*ts**(2*c-1-s)*sp.Mul(*[K(t[r],ts) for r in S if r<s])
    for s in sorted(S,reverse=True): expr=sp.residue(sp.together(expr),t[s],1)
    return sp.simplify(expr)
print('--- decoupling check')
for q in range(4):
    for S in itertools.combinations([1,2,3],q):
        pred=sp.Mul(*[-sp.Rational((-1)**lam[r-1],2**(lam[r-1]+1)) for r in [1,2,3] if r not in S])*(-1)**len(S)*(JS(lam,c,S) if S else 1)
        print(S, term(lam,c,S), pred)
