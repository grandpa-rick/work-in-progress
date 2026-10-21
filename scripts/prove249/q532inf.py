# u->oo peel of Psi with consecutive top block: Clio's single child at step 1, and the recurring head/tail tie
import sympy as sp, random
from fractions import Fraction as Fr
from q532 import Psi, pt
U=sp.Symbol('u'); V=sp.Symbol('v')
def lead_inf(kap,th,M,x):
    n=M[-1]; xx=dict(x); xx[n]=1/V          # u = 1/v
    xx={k:(sp.Rational(val.numerator,val.denominator) if isinstance(val,Fr) else val) for k,val in xx.items()}
    e=sp.cancel(sp.together(Psi(kap,th,M,xx)))
    num,den=sp.fraction(e); pn=sp.Poly(num,V); pd=sp.Poly(den,V)
    on=min(m[0] for m in pn.monoms()); od=min(m[0] for m in pd.monoms())
    cn=[c for (m,),c in zip(pn.monoms(),pn.coeffs()) if m==on][0]; cd=[c for (m,),c in zip(pd.monoms(),pd.coeffs()) if m==od][0]
    return -(on-od), cn/cd    # u-degree, leading coeff
def R(x): return {k:sp.Rational(v.numerator,v.denominator) for k,v in x.items()}
def Fth(lam): return tuple(lam[j]-2*(j+1) for j in range(len(lam)))
print("== step 1: F-block lambda (Psi twist lam_j-2j), u->oo; Clio: unique head r, child C*Psi_{lam-, const e} ==")
for lam,M in [((2,3,4),(2,3,4,5,6)),((3,4),(1,3,4,5)),((2,3),(2,3,4,5)),((1,3,4),(1,2,3,4,5,6)),((3,4,5),(1,3,4,5,6))]:
    random.seed(hash((lam,M))%1000); x=pt(max(M))
    deg,c=lead_inf(lam,Fth(lam),M,x)
    p=len(lam)-1
    while p>0 and lam[p-1]==lam[p]-1: p-=1
    b=len(lam)-p; t0=lam[p]-2*(p+1); e=t0-(b-1)
    C=sp.Integer(1)
    for m in range(b-1): C*=R(x)[lam[p]+m]**(b-1-m)
    chth=list(Fth(lam)[:p])+[e]*(b-1)
    child=Psi(lam[:-1],chth,M[:-1],R(x))
    Fm=Psi(lam[:-1],Fth(lam[:-1]),M[:-1],R(x))
    print(lam,'M',M,'deg',deg,'pred',e,'| lead/(C*child)=',sp.nsimplify(c/(C*child)),'| child twist',tuple(chth),'vs F twist of lam-',Fth(lam[:-1]))
print("== step 2+: constant-twist block node, head order e vs tail -R: tie iff e=-R ==")
for kap,th,M in [((2,3),(-2,-2),(2,3,4,5)),((2,),(-1,),(2,3,4)),((2,3),(-1,-1),(2,3,4,5)),((3,4),(-2,-2),(1,3,4,5)),((3,),(-1,),(1,3,4))]:
    random.seed(7); x=pt(max(M)); r=len(kap)
    deg,c=lead_inf(kap,th,M,x)
    head=Psi(kap[:-1],th[:-1],M[:-1],R(x)) if r>1 else sp.Integer(1)
    tail=(-1)**r*Psi(kap,tuple(t+1 for t in th),M[:-1],R(x))
    e=th[-1]
    if e==-r: pred=head+tail; lab='TIE'
    elif e>-r: pred=head; lab='head'
    else: pred=tail; lab='tail'
    print(kap,th,M,lab,'deg',deg,'pred',max(e,-r),'lead/pred=',sp.nsimplify(c/pred))
