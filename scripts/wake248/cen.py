"""Sign-coherence census for the A-residues R^A(lam,mu,nu;a,b;c), exact over Q(w) (flint).
Usage: python3 cen.py N c1,c2,...   -> writes cohN_c<c>.txt lines (one per allowed case) and a summary line.
Per case: R = Num/Den in lowest terms.
  P1 (coefficient coherence): all nonzero coefficients of Num have one sign (orientation of the
      Weyl-type factors of Den only flips a global sign, so this is well defined).
  WEYL: every irreducible factor of Den is one of w_i, w_i-1, w_i+1, w_i-w_j, w_i w_j-1 (else flagged).
  P1* (MAIN, Weyl-stripped coherence): Num = monomial * prod(C_n-root factors)^e * P with P coprime to them;
      coherent iff all coefficients of P have one sign.
  P3 (positive-orthant evaluation coherence): sign of Num at 40 random rational points with all w>0.
  sgn: the common sign (+1/-1) when P1 holds."""
import sys,itertools,random,time
from fractions import Fraction
from eng import *
N=int(sys.argv[1]); cs=list(map(int,sys.argv[2].split(',')))
C=setup(N); Ss=subsets(N); pairs=list(itertools.combinations(range(1,N+1),2))
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
R_=random.Random(7)
pts=[[Fraction(R_.randint(1,40),R_.randint(1,40)) for _ in range(N)] for _ in range(40)]
def analyse(r):
    if r.n==0: return dict(zero=True)
    co=[int(v*1) if False else v for v in r.n.to_dict().values()]
    sg=set((v>0)-(v<0) for v in co)
    p1=len(sg)==1
    facs=[f for f,e in r.d.factor()[1]]
    weyl=all(weyl_factor(f) for f in facs)
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
import flint
cases=[]
for i,l in enumerate(Ss):
    for m in Ss[i:]:
        for k in Ss:
            for a,b in pairs:
                if EA(N,a,b,l,m,k): cases.append((a,b,l,m,k))
syss={p:Sys(N,*p) for p in pairs}
for c in cs:
    t0=time.time(); out=open('coh%d_c%d.txt'%(N,c),'w'); cnt=dict(cohS=0,incS=0,coh=0,inc=0,zero=0,nonweyl=0,p3inc=0,plus=0,minus=0)
    for (a,b,l,m,k) in cases:
        r=A(syss[a,b],a,b,l,m,k,c); d=analyse(r)
        out.write(repr(((a,b,l,m,k),d))+'\n')
        if d['zero']: cnt['zero']+=1; continue
        cnt['coh' if d['p1'] else 'inc']+=1
        cnt['cohS' if d['p1s'] else 'incS']+=1
        if not d['weyl']: cnt['nonweyl']+=1
        if not d['p3']: cnt['p3inc']+=1
        if d['p1']: cnt['plus' if d['sgn']>0 else 'minus']+=1
    line='N %d c %d allowed %d %s %.0fs'%(N,c,len(cases),cnt,time.time()-t0)
    out.write('# '+line+'\n'); out.close(); print(line,flush=True)
