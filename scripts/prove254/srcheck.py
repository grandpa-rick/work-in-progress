# PROVE 254: independent check of SR_LG / SR_OG: Demazure chi at NEGATIVE j vs brute-force strict counts on P_kappa.
import sys, itertools
from functools import lru_cache
sys.path.insert(0,'../prove249d'); sys.path.insert(0,'../wake253')
from dem import chi as chiC
import og
def boxes(k): return sorted({(i,i+j) for i,li in enumerate(k) for j in range(li)}, key=lambda b:(b[0]+b[1],b))
def count(k, nv, diag_even):
    P=boxes(k); idx={b:t for t,b in enumerate(P)}
    preds=[[idx[q] for q in ((a-1,b),(a,b-1)) if q in idx] for (a,b) in P]   # covers suffice (product order on shifted diag)
    f=[0]*len(P); tot=0
    def rec(t):
        nonlocal tot
        if t==len(P): tot+=1; return
        lo=max([f[s] for s in preds[t]],default=0)+1
        for v in range(lo,nv+1):
            if diag_even and P[t][0]==P[t][1] and v%2: continue
            f[t]=v; rec(t+1)
    rec(0); return tot
n=int(sys.argv[1]); bad=0; tot=0
for L in range(n+1):
  for c in itertools.combinations(range(1,n+1),L):
    k=tuple(sorted(c,reverse=True))
    for j in range(1,5):
        hl=chiC(set(k),n,-j)[0]; hl_pred=(-1)**sum(k)*count(k,2*j-1,True)
        ho=og.chi_dem(k,n,-j); ho_pred=(-1)**sum(k)*count(k,j-1,False)
        tot+=1
        if hl!=hl_pred or ho!=ho_pred: bad+=1; print('FAIL',k,j,hl,hl_pred,ho,ho_pred)
print('n',n,'checked',tot,'bad',bad)
