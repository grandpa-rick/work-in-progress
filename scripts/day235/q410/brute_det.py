#!/usr/bin/env python3
"""Brute-force t-level check of warnaar_calpha.py: expand Warnaar (6.1) as an honest Laurent polynomial
in n variables (sympy) and read off [x^alpha] (x1..xn)^k G. Small cases only."""
import itertools, sympy as sp
from warnaar_calpha import calpha, t
def run(k, l, n, alphas):
    M=k+l+1; N=2*M; xs=sp.symbols(f'x1:{n+1}')
    vals=[v for x in xs for v in (x,1/x)]
    u=sp.symbols('u'); gen=sp.expand(sp.prod([1+u*v for v in vals]))
    E=lambda r: gen.coeff(u,r) if 0<=r<=2*n else 0
    G=0
    yr=range(-(2*n//N)-2,(2*n//N)+3)
    for y in itertools.product(yr,repeat=k):
        mat=sp.Matrix(k,k,lambda i,j: t**(M*y[i]**2-(j+1)*y[i])*(E(n-(i+1)+(j+1)-N*y[i])-E(n+(i+1)+(j+1)-N*y[i])))
        G+=mat.det()
    F=sp.expand(G*sp.prod(xs)**k)
    P=sp.Poly(F,*xs)
    for a in alphas:
        b=sp.expand(P.coeff_monomial(sp.prod([x**e for x,e in zip(xs,a)])))
        c=calpha(k,l,list(a))
        print(f"(k,l,n)=({k},{l},{n}) alpha={a}: brute={b}  engine={c}  equal={sp.expand(b-c)==0}")
run(1,0,4,[(1,1,1,1),(2,0,1,1),(0,0,2,2),(2,2,0,0)])
run(1,1,4,[(1,1,1,1),(2,1,0,1)])
run(2,0,3,[(1,1,1),(2,1,0),(4,2,0)])
