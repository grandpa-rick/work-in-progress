# Theorem S (symbol theorem) exact test.  u = w_n symbolic, other w random rationals.
import sympy as sp, itertools, random
random.seed(249)
u,z=sp.symbols('u z')
def setup(N,n):
    w={k:sp.Rational(random.randint(2,60),random.randint(2,60))+k for k in range(1,N+1)}
    w[n]=u; return w
def S0(M,H,w):
    x=[w[h] for h in H]; e=sp.Integer(1)
    for a in x: e*=(1-a**-2)
    for i in range(len(x)):
        for j in range(i+1,len(x)): e*=(1-1/(x[i]*x[j]))
    for t in M:
        if t in H: continue
        for a in x: e*=(a-1/w[t])/(a-w[t])
    return e
def Phi(M,r,f,w):   # f: function of list of r values (symmetric)
    s=0
    for H in itertools.combinations(M,r):
        s+=S0(M,H,w)*f([w[h] for h in H])
    return s
def rho(k,y,w):
    p=sp.Integer(1)
    for i in range(1,k): p*=(1-w[i]/y)
    return p
def Aalt(kap,th,w):
    r=len(kap); X=sp.symbols('X0:%d'%r)
    M=sp.Matrix(r,r,lambda i,j: X[i]**(j+th[j])*rho(kap[j],X[i],w))
    V=sp.Integer(1)
    for i in range(r):
        for j in range(i+1,r): V*=(X[j]-X[i])
    A=sp.cancel(sp.together(M.det()/V))
    return A,X
def lead(expr):
    expr=sp.cancel(sp.together(expr)); nu,de=sp.fraction(expr)
    pn=sp.Poly(nu,u); pd=sp.Poly(de,u)
    on=min(m[0] for m in pn.monoms()); od=min(m[0] for m in pd.monoms())
    return on-od, sp.nsimplify(pn.coeff_monomial(u**on)/pd.coeff_monomial(u**od))
def test(kap,th,M):
    N=max(M); n=N; w=setup(N,n); r=len(kap); Mm=tuple(m for m in M if m!=n)
    Psi=Phi(M,r,lambda xs: Aalt(kap,th,w)[0].subs(dict(zip(Aalt(kap,th,w)[1],xs))),w) if False else None
    A,X=Aalt(kap,th,w)
    fA=lambda xs: A.subs(dict(zip(X,xs)))
    o_true,c_true=lead(Phi(M,r,fA,w))
    # symbol: z-order of A(x',z) with x'=X[:-1] symbols, z=X[-1]
    Az=sp.Poly(sp.expand(A*X[-1]**40),X[-1])   # shift to make polynomial
    o=min(m[0] for m in Az.monoms())-40
    ao=sp.expand(Az.coeff_monomial(X[-1]**(o+40)))
    eM=sp.Mul(*[w[m] for m in Mm])
    def head(oo):
        a=sp.expand(Az.coeff_monomial(X[-1]**(oo+40))) if oo+40>=0 else 0
        f=lambda xs: (sp.Mul(*xs)*a).subs(dict(zip(X[:-1],xs)))
        return -eM**-2*(-1)**(r-1)*Phi(Mm,r-1,f,w) if r>1 else -eM**-2*(-1)**(r-1)*a
    tail=lambda: (-1)**r*Phi(Mm,r,lambda xs:(A/sp.Mul(*X)).subs(dict(zip(X,xs))),w)
    if o<1: o_pred,c_pred=o-r-1,head(o)
    elif o>1: o_pred,c_pred=-r,tail()
    else: o_pred,c_pred=-r,tail()+head(1)
    ok=(o_pred==o_true) and sp.simplify(c_pred-c_true)==0
    print(kap,th,M,'z-order o=',o,'pred order',o_pred,'true',o_true,'coeff match',ok)
    return ok
for kap,th,M in [((1,3),(0,1),(1,2,3,4)),((1,3),(2,2),(1,2,3,4)),((1,3),(1,2),(1,2,3,4)),((1,2),(1,2),(1,2,3)),
                 ((2,4),(0,1),(1,2,4,5)),((1,2,4),(0,0,1),(1,2,3,4,5)),((1,3,4),(0,1,1),(1,3,4,5)),((2,4),(3,3),(2,3,4,5)),((1,3),(3,1),(1,2,3,4))]:
    test(kap,th,M)
