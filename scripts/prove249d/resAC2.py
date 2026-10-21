# general pair: find nu~ subset [N-2] with Res rho_{nu,L''+ab} = K w_{L''}^{-2} rho^{(N-2)}_{nu~,L''} (K const in L'')
import sys, random
from fractions import Fraction as Fr
sys.path.insert(0,'../prove249c')
from eng import ell, leq, subsets
N=int(sys.argv[1]); random.seed(4)
base={i:Fr(random.randint(2,40),random.randint(2,40)) for i in range(1,N+1)}
def R(nu,L,a,b):
    out=[]
    for dl in (Fr(1,10**25),Fr(2,10**25)):
        w=dict(base); w[b]=1/(w[a]*(1+dl)); h=1-1/(w[a]*w[b])
        out.append(ell(nu,L,N,w)*h)
    return 2*out[0]-out[1]
tot=found=0; Ks={}
for a in range(1,N+1):
  for b in range(a+1,N+1):
    rest=[i for i in range(1,N+1) if i not in (a,b)]
    wsm={i+1:base[x] for i,x in enumerate(rest)}
    for nu in subsets(N):
        Ls=[Lpp for Lpp in subsets(N-2) if leq(sorted([rest[i-1] for i in Lpp]+[a,b]),nu)]
        if not Ls: continue
        tot+=1
        vals={tuple(Lpp):R(nu,sorted([rest[i-1] for i in Lpp]+[a,b]),a,b) for Lpp in Ls}
        hit=None
        for nt in subsets(N-2):
            ok=True; K=None
            for Lpp in subsets(N-2):
                v=vals.get(tuple(Lpp),Fr(0))
                s=ell(nt,Lpp,N-2,wsm) if leq(Lpp,nt) else Fr(0)
                p=Fr(1)
                for i in Lpp: p*=wsm[i]**2
                if s==0:
                    if abs(v)>Fr(1,10**15): ok=False;break
                    continue
                k=v*p/s
                if K is None: K=k
                elif abs(k-K)>abs(K)*Fr(1,10**12): ok=False;break
            if ok and K is not None: hit=(tuple(nt),float(K)); break
        if hit: found+=1; Ks.setdefault((a,b),set()).add(round(hit[1],8))
        else: print('no nu~ for',(a,b),nu, 'support',Ls)
print('N',N,'instances',tot,'with clean nu~',found); print({k:v for k,v in Ks.items()})
