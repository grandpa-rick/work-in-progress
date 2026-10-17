"""N=4 control: all 127 jump>=2 triples (lam<=mu in order), every pair (a,b) in Clio's predicted set (Conj divisor);
c_1 residue on H_ab via day243/q488c/resnum.residue (proved first-order residue formula, numeric w, 2 seeds).
Stratum (Clio 10-17 reply 4.4): nonvanH instances (M=kappa\\{a,b}, I, a, b) with I in {lam, mu}; in-stratum iff max I <= b.
Pair-level 'in' iff max(lam u mu) <= b (both labels in stratum); also report mu-only / lam-only."""
import sys, os, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../day243/q488c'))
sys.argv=[sys.argv[0],'1']
from engine import subsets, geq
from coverage import cover
from resnum import residue
N=4; S=subsets(N)
def rm(M,a,b): return tuple(x for x in M if x not in (a,b))
def pairs(l,m,nu):
    P=[]
    for k in S:
        if not geq(nu,k): continue
        for a,b in itertools.combinations(k,2):
            r=rm(k,a,b)
            if geq(r,l) and geq(r,m) and (a,b) not in P: P.append((a,b))
    return P
tot,ct,cs,unc=cover(N)
T=[]
for i,l in enumerate(S):
    for m in S[i:]:
        for nu in S:
            if geq(nu,l) and geq(nu,m) and len(nu)-max(len(l),len(m))>=2: T.append((l,m,nu))
print('triples',len(T))
mx=lambda X: max(X) if X else 0
stat={}  # (in_stratum, c1pole) -> count over (triple,pair)
trip={}
for (l,m,nu) in T:
    P=pairs(l,m,nu); hits=[]
    for a,b in P:
        h=residue(l,m,nu,a,b,N,seed=1)!=0 or residue(l,m,nu,a,b,N,seed=7)!=0
        ins=mx(l+m)<=b
        stat[ins,h]=stat.get((ins,h),0)+1
        hits.append(((a,b),'in' if ins else 'out','c1pole' if h else 'noc1pole'))
    succ=any(x[2]=='c1pole' for x in hits)
    anyin=any(x[1]=='in' for x in hits)
    trip[anyin,succ]=trip.get((anyin,succ),0)+1
    print(l,m,nu,'maxlam',mx(l),'maxmu',mx(m),hits,'SUCCESS' if succ else 'C-FAIL',flush=True)
print('PAIR-LEVEL (in_stratum,c1pole):',stat)
print('TRIPLE-LEVEL (some predicted pair in-stratum, c-success):',trip)
