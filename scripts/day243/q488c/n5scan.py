import itertools, sys, time
from engine import subsets, geq
from coverage import cover
from resnum import residue
N=int(sys.argv[1]); t0=time.time()
tot,ct,cs,unc=cover(N)
S=subsets(N); fails=[]
def rm(M,a,b): return tuple(x for x in M if x not in (a,b))
for (l,m,nu) in unc:
    pairs=[]
    for k in S:
        if not(geq(nu,k)): continue
        for a,b in itertools.combinations(k,2):
            r=rm(k,a,b)
            if geq(r,l) and geq(r,m) and (a,b) not in pairs: pairs.append((a,b))
    hit=None
    for a,b in pairs:
        if residue(l,m,nu,a,b,N,seed=3)!=0: hit=(a,b); break
    print(l,m,nu,'q1-pole at',hit if hit else 'NONE (c fails)', len(pairs),'pairs tried',flush=True)
    if not hit: fails.append((l,m,nu))
print('uncovered',len(unc),'c-fails',len(fails),'time',time.time()-t0)
