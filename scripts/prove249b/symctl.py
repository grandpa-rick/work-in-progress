exec(open('sym.py').read().split('for kap,th,M in')[0])
def ctl(kap,th,M,mode):
    N=max(M); n=N; w=setup(N,n); r=len(kap); Mm=tuple(m for m in M if m!=n)
    A,X=Aalt(kap,th,w); fA=lambda xs: A.subs(dict(zip(X,xs)))
    o_true,c_true=lead(Phi(M,r,fA,w))
    Az=sp.Poly(sp.expand(A*X[-1]**40),X[-1]); eM=sp.Mul(*[w[m] for m in Mm])
    a=sp.expand(Az.coeff_monomial(X[-1]**41))
    hd=-eM**-2*(-1)**(r-1)*Phi(Mm,r-1,lambda xs:(sp.Mul(*xs)*a).subs(dict(zip(X[:-1],xs))),w)
    tl=(-1)**r*Phi(Mm,r,lambda xs:(A/sp.Mul(*X)).subs(dict(zip(X,xs))),w)
    pred={'tail only':tl,'head only':hd,'head sign flipped':tl-hd}[mode]
    print(kap,th,mode,'breaks:',sp.simplify(pred-c_true)!=0)
for kap,th,M in [((1,3),(2,2),(1,2,3,4)),((1,3),(1,2),(1,2,3,4)),((2,4),(3,3),(2,3,4,5))]:
    for mode in ['tail only','head only','head sign flipped']: ctl(kap,th,M,mode)
