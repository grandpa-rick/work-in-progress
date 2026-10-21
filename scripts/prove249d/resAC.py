# A->C equivariant probe: Res_{H_ab} rho_{nu,L''+ab} vs rho^{(N-2)}_{nu--, L''}
import sys, random, itertools
from fractions import Fraction as Fr
sys.path.insert(0,'../prove249c')
from eng import ell, leq, subsets
N=int(sys.argv[1]); a,b=int(sys.argv[2]),int(sys.argv[3]); random.seed(4)
base={i:Fr(random.randint(2,40),random.randint(2,40)) for i in range(1,N+1)}
def resrho(nu,L,dl):
    w=dict(base); w[b]=1/(w[a]*(1+dl))
    h=1-1/(w[a]*w[b])
    return ell(nu,L,N,w)*h
def R(nu,L):
    d1,d2=Fr(1,10**25),Fr(2,10**25)
    v1,v2=resrho(nu,L,d1),resrho(nu,L,d2)
    return (2*v1-v2)   # Richardson O(d^2)
rest=[i for i in range(1,N+1) if i not in (a,b)]
ren={x:i+1 for i,x in enumerate(rest)}
wsm={ren[x]:base[x] for x in rest}
for nu in subsets(N):
    if not (a in nu and b in nu) and True: pass
    if len(nu)<2: continue
    numm=sorted(nu)[:-2]
    out=[]
    for Lpp in subsets(N-2):
        L=sorted([rest[i-1] for i in Lpp]+[a,b])
        if not leq(L,nu): continue
        r=R(nu,L)
        s=ell(numm,Lpp,N-2,wsm) if leq(Lpp,numm) else Fr(0)
        out.append((tuple(Lpp), float(r), float(s), float(r/s) if s!=0 else None))
    print(nu,'->',tuple(numm), out)
print('--- n(L'') * w_{L''}^2 :')
nu=list(range(1,N+1)); numm=sorted(nu)[:-2]
for Lpp in subsets(N-2):
    L=sorted([rest[i-1] for i in Lpp]+[a,b])
    r=R(nu,L); s=ell(numm,Lpp,N-2,wsm)
    p=Fr(1)
    for i in Lpp: p*=wsm[i]**2
    print(tuple(Lpp), float(r/s*p))
