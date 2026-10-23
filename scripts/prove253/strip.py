# c_rho != 0 for rho* = kappa minus top t antidiagonals (t=1 r odd, t=2 r even), via BR's h-sum formula (brpieri.C), and sign rule of all c.
import sys, itertools
sys.path.insert(0,'../wake252/caseB'); sys.path.insert(0,'../wake251/orderpoly')
from brpieri import pieri_dim, boxes
from orderpoly import shifted, boxless, longest_chain
def sp(n): return [tuple(sorted(c,reverse=True)) for L in range(n+1) for c in itertools.combinations(range(1,n+1),L)]
def r(k): return longest_chain(shifted(list(k)),boxless) if k else 0
def frombox(B):
    rows={}
    for (a,b) in B: rows[a]=rows.get(a,0)+1
    return tuple(rows[a] for a in sorted(rows))
if __name__=="__main__":
 for n in range(1,int(sys.argv[1])+1):
     ok=sg=tot=rf=0; bad=[]
     for k in sp(n):
         if not k: continue
         tot+=1; R=r(k); B=boxes(k)
         rf += R==max(k[i]+2*i for i in range(len(k)))
         t=1 if R%2 else 2
         rs=frombox({x for x in B if x[0]+x[1] < R-t})
         P=pieri_dim(k,n)
         if P.get(rs,0)*(-1)**(sum(k)-sum(rs)-1) > 0: ok+=1
         else: bad.append((k,rs,P.get(rs,0)))
         sg += all(v*(-1)**(sum(k)-sum(p)-1) > 0 for p,v in P.items())
     print(f'm={n}: r=max(k_i+2(i-1)) {rf}/{tot}; c(rho*) nonzero with sign (-1)^(|th|-1): {ok}/{tot} {bad[:4]}; all c sign-alternating: {sg}/{tot}')
