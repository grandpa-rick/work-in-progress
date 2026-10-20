import random,flint
from fractions import Fraction
def weyl_factor(f):
    d=f.degrees(); coeffs=f.to_dict()
    tot=sum(d)
    if tot==0: return True
    mons=sorted(coeffs.items())
    if len(mons)==1 and tot==1: return True          # w_i
    if len(mons)==2:
        (e1,c1),(e2,c2)=mons
        s1,s2=sum(e1),sum(e2)
        if abs(c1)==1 and abs(c2)==1:
            if {s1,s2}=={0,1}: return True               # w_i +- 1
            if s1==s2==1 and c1*c2<0: return True         # w_i - w_j
            if {s1,s2}=={0,2} and max(max(e1),max(e2))==1 and c1*c2<0: return True  # w_i w_j - 1
    return False
def analyse(r,N):
    if r.n==0: return dict(zero=True)
    co=[int(v*1) if False else v for v in r.n.to_dict().values()]
    sg=set((v>0)-(v<0) for v in co)
    p1=len(sg)==1
    facs=[f for f,e in r.d.factor()[1]]
    weyl=all(weyl_factor(f) for f in facs)
    R_=random.Random(7); pts=[[Fraction(R_.randint(1,40),R_.randint(1,40)) for _ in range(N)] for _ in range(40)]
    nf=r.n
    ev=set()
    for p in pts:
        v=nf(*[flint.fmpq(x.numerator,x.denominator) for x in p])
        ev.add((v>0)-(v<0))
    ev.discard(0)
    # P1*: strip every C_n-root factor (w_i, w_i+-1, w_i-w_j, w_iw_j-1) from Num; test the residual P
    cst,nf_=r.n.factor(); P=r.n.context().from_dict({(0,)*N:1})*cst; nres=0
    for f,e in nf_:
        if not weyl_factor(f): P=P*f**e; nres+=1
    sgP=set((v>0)-(v<0) for v in P.to_dict().values())
    return dict(zero=False,p1=p1,sgn=(list(sg)[0] if p1 else 0),weyl=weyl,p3=(len(ev)==1),nterms=len(co),
                p1s=(len(sgP)==1),Pterms=len(P.to_dict()),nres=nres,P=str(P) if len(sgP)>1 else '')
