# s_{a+1/2}(k) = sum_{L=0}^{k1} c_L N_{k/(L)}(a), c_0=1, c_L=-(-1/2)^L ; N_{k/(L)}(a)=#strict maps of P_k minus first L boxes of row 0 -> [2a-1], diag even
import sys, itertools
from fractions import Fraction as Fr
sys.argv=[sys.argv[0]]+sys.argv[1:]
import half
def cntskew(k,L,nv):
    P=sorted({(i,i+j) for i,li in enumerate(k) for j in range(li)}-{(0,j) for j in range(L)}, key=lambda b:(b[0]+b[1],b))
    idx={b:t for t,b in enumerate(P)}
    preds=[[idx[q] for q in ((a-1,b),(a,b-1)) if q in idx] for (a,b) in P]
    f=[0]*len(P); tot=0
    def rec(t):
        nonlocal tot
        if t==len(P): tot+=1; return
        lo=max([f[s] for s in preds[t]],default=0)+1
        for v in range(lo,nv+1):
            if P[t][0]==P[t][1] and v%2: continue
            f[t]=v; rec(t+1)
    rec(0); return tot
N=int(sys.argv[1]); A=int(sys.argv[2])
bad=0; tot=0; ineq_fail=[]
for k in half.KS:
    if not k: continue
    for a in range(1,A+1):
        Ns=[cntskew(k,L,2*a-1) for L in range(k[0]+1)]
        v=sum((Fr(1) if L==0 else -Fr(-1,2)**L)*Ns[L] for L in range(k[0]+1))
        tot+=1
        if v!=half.s(k,2*a+1): bad+=1; print('FORMULA FAIL',k,a,v,half.s(k,2*a+1))
        if 2*a-1>=k[0]:
            for L in range(1,k[0],2):
                if Ns[L+1]>2*Ns[L]: ineq_fail.append((k,a,L,Ns))
print('formula checked',tot,'bad',bad)
print('pair inequality N_{L+1}<=2N_L (L odd, 2a-1>=k1) failures:',len(ineq_fail), ineq_fail[:8])
