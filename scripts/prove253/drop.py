import sys, itertools
sys.path.insert(0,'../wake251/orderpoly')
from orderpoly import shifted, boxless, longest_chain
def sp(n): return [tuple(sorted(c,reverse=True)) for L in range(n+1) for c in itertools.combinations(range(1,n+1),L)]
def boxes(l): return {(i,i+j) for i,li in enumerate(l) for j in range(li)}
def r(k): return longest_chain(shifted(list(k)),boxless) if k else 0
for n in range(1,int(sys.argv[1])+1):
    S=sp(n); B={k:boxes(k) for k in S}; R={k:r(k) for k in S}
    cL=cO=bL=bO=0; mL=mO=0
    for k in S:
        for p in S:
            if p==k or not B[p]<=B[k]: continue
            th=B[k]-B[p]
            rim=not any(a<c and b<d for (a,b) in th for (c,d) in th)
            rows=max([sum(1 for x in th if x[0]==i) for i in range(n)]); cols=max([sum(1 for x in th if x[1]==j) for j in range(n)])
            if rim and rows<=2 and cols<=2:
                cL+=1; bL+= R[p]<R[k]-2; mL=max(mL,R[k]-R[p])
            if rim and rows<=1 and cols<=1:
                cO+=1; bO+= R[p]<R[k]-1; mO=max(mO,R[k]-R[p])
    print(f'm={n}: LG-admissible pairs {cL}, violations of r(rho)>=r(k)-2: {bL}, max drop {mL} | OG-admissible {cO}, violations of r(rho)>=r(k)-1: {bO}, max drop {mO}')
