import sys, itertools, sympy as sp
from span import G
def hpoly(n,xs):
    if n<0: return 0
    return sum(sp.Mul(*c) for c in itertools.combinations_with_replacement(xs,n)) if n>0 else 1
ok=0;tot=0
for k in range(1,6):
  for a in range(1,5):
    y=sp.symbols('y1:%d'%(k+1)); x=[1/v for v in y]
    D=sp.expand(sum(c*sp.Mul(*[y[j]**(-kk[j]) for j in range(k)]) for kk,c in G(k,a).items()))
    ew=sp.Mul(*[v**a for v in y]); F=sp.Mul(*[v**(a-1)*(1+v) for v in y])
    pred=ew-(-1)**(k-1)*F*sum(hpoly(2*s-k+1,x) for s in range(a))
    tot+=1; ok+= sp.simplify(D-pred)==0
print('G_k(a) formula',ok,'/',tot)
