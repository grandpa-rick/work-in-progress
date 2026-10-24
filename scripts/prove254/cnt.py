# PROVE 254: independent check of SR_LG / SR_OG: Demazure chi at NEGATIVE j vs brute-force strict counts on P_kappa.
import sys, itertools
from functools import lru_cache
sys.path.insert(0,'../prove249d'); sys.path.insert(0,'../wake253')
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
