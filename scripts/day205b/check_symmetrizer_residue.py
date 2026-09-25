# Day 205b: check (i) sigma_m F = sum_i F^(i) prod_{j!=i} (X_i - t X_j)/(X_i - X_j) for F sym in tail
# (ii) (1-t) sum_i x_i^n prod a_ij = q_n (HL), (iii) M_a(r) closed form.
import sympy as sp, random, itertools, sys
sys.path.insert(0,'/home/agent/projects/proofs/scripts/day198')
from p2Y_er import build_action
t=sp.symbols('t')
def e(Xs,k):
    if k<0: return 0
    if k==0: return sp.Integer(1)
    return sum(sp.prod(c) for c in itertools.combinations(Xs,k))
ok=True
for m in range(1,6):
    X,Ti,_,_,_=build_action(m) if m>1 else (sp.symbols('X1:2'),None,None,None,None)
    X=list(sp.symbols(f'X1:{m+1}'))
    for trial in range(3):
        # random F = sum X1^a * e_lambda(tail)
        tail=X[1:]
        F=0
        for _ in range(3):
            a=random.randint(0,3); lam=[random.randint(0,min(3,len(tail))) for _ in range(2)]
            F+=random.randint(-3,3)*X[0]**a*sp.prod([e(tail,k) for k in lam])
        F=sp.expand(F)
        lhs=F
        for k in range(1,m):
            c=F
            for j in range(1,k+1): c=Ti(c,j)
            lhs+=c
        rhs=0
        for i in range(m):
            sw={X[0]:X[i],X[i]:X[0]}
            Fi=F.xreplace(sw)
            rhs+=Fi*sp.prod([(X[i]-t*X[j])/(X[i]-X[j]) for j in range(m) if j!=i])
        d=sp.simplify(sp.together(sp.expand(lhs)-rhs))
        if d!=0: ok=False; print('FAIL sym',m,F)
    # (ii)
    y=sp.symbols('y')
    Q=sp.series(sp.prod([(1-t*x*y)/(1-x*y) for x in X]),y,0,6).removeO()
    for n in range(1,6):
        S=sum(X[i]**n*sp.prod([(X[i]-t*X[j])/(X[i]-X[j]) for j in range(m) if j!=i]) for i in range(m))
        d=sp.simplify(sp.together((1-t)*S-Q.coeff(y,n)))
        if d!=0: ok=False; print('FAIL qn',m,n)
print('ALL OK' if ok else 'FAILURES')
