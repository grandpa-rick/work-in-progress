import sympy as sp
from k4_analyze import load, qint, q, t, u
data=load()
rs=sorted({r for r,m in data})
for r in rs:
    tau=sp.cancel(data[(r,min(m for rr,m in data if rr==r))])
    n,d=sp.fraction(tau)
    print(f'\nr={r}: denom={d}')
    c,fl=sp.factor_list(n)
    for fac,e in fl:
        print(f'   deg(q,t)={sp.Poly(fac,q,t).degree_list()} ^{e}:', fac if len(str(fac))<70 else str(fac)[:70]+'...')
    Q=sp.cancel(n*qint(4)/((q**4-1)*qint(r+4)))
    qn,qd=sp.fraction(Q)
    print('   quotient denom (should be 1 or Phi_4/Phi_2 leftover):', sp.factor(qd))
