import sys, itertools, collections
sys.path.insert(0,'.')
import brpieri as B
from support import rp, bx
for n in range(2, int(sys.argv[1])+1):
    strict=[tuple(sorted(c,reverse=True)) for L in range(n+1) for c in itertools.combinations(range(1,n+1),L)]
    stair=bx(tuple(range(n,0,-1)))
    R=lambda s:{(n-1-c,n-1-r) for (r,c) in s}
    refl=all(bx(B.dual(k,n))==R(stair-bx(k)) for k in strict)
    bad=0; drops=collections.Counter(); tot=0
    for k in strict:
        P=B.pieri_dim(k,n)
        for rho in strict:
            th=bx(k)-bx(rho)
            if not (len(rho)<=len(k) and all(rho[i]<=k[i] for i in range(len(rho)))) or rho==k: 
                if rho in P: bad+=1
                continue
            rim=not any(r2>r1 and c2>c1 for (r1,c1) in th for (r2,c2) in th)
            rc=max(collections.Counter(a for a,_ in th).values()) <=2 and max(collections.Counter(b for _,b in th).values())<=2
            if (rho in P) and not (rim and rc): bad+=1
        if len(k)>=2:
            tot+=1; drops[min(rp(r)-rp(k) for r in P)]+=1
    print(n,'reflection',refl,'support-rule mismatches',bad,'min drop dist',sorted(drops.items()),'of',tot)
