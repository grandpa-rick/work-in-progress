exec(open('sym.py').read().split('for kap,th,M in')[0])
def R2(kap,th,M):
    N=max(M); n=N; w=setup(N,n); r=len(kap); Mm=tuple(m for m in M if m!=n)
    A,X=Aalt(kap,th,w); o_true,c_true=lead(Phi(M,r,lambda xs:A.subs(dict(zip(X,xs))),w))
    b=sp.expand(sp.cancel(A/sp.Mul(*X)))
    assert sp.Poly(b,*X).monoms() and min(min(m) for m in sp.Poly(b,*X).monoms())>=0
    b0=b.subs(X[-1],0)           # function of X[:-1]
    def Lag(F,xs):
        s=0
        for t in range(r):
            oth=[xs[j] for j in range(r) if j!=t]
            pr=sp.Mul(*[xs[j]/(xs[j]-xs[t]) for j in range(r) if j!=t])
            s+=F.subs(dict(zip(X[:-1],oth)),simultaneous=True)*pr
        return sp.cancel(s)
    LB=Lag(b0,list(X)); B=sp.expand(b-LB)
    Bp=sp.cancel(B/sp.Mul(*X)); divisible=sp.fraction(Bp)[1].free_symbols.isdisjoint(set(X))
    pred=(-1)**r*(Phi(Mm,r-1,lambda xs:b0.subs(dict(zip(X[:-1],xs))),w)+Phi(Mm,r,lambda xs:B.subs(dict(zip(X,xs))),w))
    print(kap,th,'order',o_true,'R2 match',sp.simplify(pred-c_true)==0,'B=',sp.factor(B),'e|B',divisible,'b0=',b0)
for kap,th,M in [((1,3),(2,2),(1,2,3,4)),((1,3),(1,2),(1,2,3,4)),((1,2),(1,2),(1,2,3)),((2,4),(3,3),(2,3,4,5)),((1,2,4),(1,1,2),(1,2,3,4,5)),((2,3),(2,2),(2,3,4))]:
    R2(kap,th,M)
