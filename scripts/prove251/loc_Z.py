import sys, itertools, sympy as sp
from span import G
def chi_loc(k,a,conv):
    y=sp.symbols('y1:%d'%(k+1)); E=lambda i: 1/y[i]   # e^{eps_i}
    tot=0
    for i in range(k):
        for top in ('F','Fi'):
            # weights of U at point (exponent vectors as lists of eps-coeffs)
            Uw=[[-1 if j==t else 0 for t in range(k)] for j in range(k)]
            if top=='Fi': Uw[i]=[1 if t==i else 0 for t in range(k)]
            mono=lambda w: sp.Mul(*[ (1/y[t])**w[t] for t in range(k)])  # e^{w}
            chU=[mono(w) for w in Uw]
            detU=sp.Mul(*chU)
            fib=sp.Mul(*[(1+c) for c in chU])*detU**a
            tw=[[ (1 if t==j else 0)-(1 if t==i else 0) for t in range(k)] for j in range(k) if j!=i]
            tw.append([ (2 if t==i else 0)*(1 if top=='F' else -1) for t in range(k)])
            den=sp.Mul(*[(1-mono(w)**(-1 if conv==0 else 1)) for w in tw])
            tot+=fib/den
    return sp.factor(sp.simplify(tot)), y
def dict2poly(f,y):
    return sp.expand(sum(c*sp.Mul(*[y[j]**(-kk[j]) for j in range(len(kk))]) for kk,c in f.items()))
for k,a in [(1,1),(2,1),(2,2),(3,2)]:
    for conv in (0,1):
        L,y=chi_loc(k,a,conv); D=dict2poly(G(k,a),y)
        print(k,a,'conv',conv,'match',sp.simplify(L-D)==0)
