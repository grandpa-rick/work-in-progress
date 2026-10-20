from psi import *
import random
# numeric check of guess: Psi_[r],0(M) + prod_M w^-2 Psi_[r-1],2(M) == Psi_[r-1],0(M) ?
def num(expr,vals): return sp.nsimplify(expr.subs(vals))
for r,M in [(2,(1,2)),(2,(1,2,3)),(2,(2,3,4)),(3,(1,2,3)),(3,(1,2,3,4))]:
    A=Psi(tuple(range(1,r+1)),(0,)*r,M); B=Psi(tuple(range(1,r)),(2,)*(r-1),M); C=Psi(tuple(range(1,r)),(0,)*(r-1),M)
    pw=sp.Mul(*[w(m)**-2 for m in M])
    print(r,M,sp.factor(A+pw*B-C), '|', sp.factor(A+pw*B))
