# PROVE 254 stretch: half-integer theorem.  s_{a+1/2}(k) = sum_g Phi(u(g)),  g: rows>=1 of P_k -> [2a-1] strict, diag even;
# Phi = E0 + sum_{i=0}^{n-2} 2^{-i-1} G_{1+i}(M-1-i, u-1-i) + 2^{-n}  if M>=n-1 ;  = (-1)^n 2^{-n} if M<=n-1  (n=k1, M=2a-1)
# Checked against Demazure (type C) interpolated to x=-a-1/2: independent of (*), BR, Thm A.
import sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0,'../prove249d')
from dem import chi as chiC
def fills(n,M,u,L,even0):
    cnt=0
    for vs in itertools.combinations(range(1,M+1),n-L) if M>=0 else []:
        ok=True
        for i,v in enumerate(vs):
            b=L+i
            if 1<=b<=len(u) and v>=u[b-1]: ok=False;break
            if b==0 and even0 and v%2: ok=False;break
        if ok: cnt+=1
    return cnt
def Phi_def(n,M,u): return fills(n,M,u,0,True)-sum(Fr(-1,2)**L*fills(n,M,u,L,False) for L in range(1,n+1))
def Phi_pos(n,M,u):
    assert M>=n-1
    return fills(n,M,u,0,True)+sum(Fr(1,2**(i+1))*fills(n,M-1-i,tuple(x-1-i for x in u),1+i,False) for i in range(n-1))+Fr(1,2**n)
def rows_maps(k,M):   # maps of rows>=1 of P_k (= P_{k'}) -> [M], diag even; yield bounds u_b=g(1,b), b=1..k2
    P=sorted({(i,i+j) for i,li in enumerate(k) for j in range(li) if i>=1}, key=lambda b:(b[0]+b[1],b))
    idx={b:t for t,b in enumerate(P)}
    preds=[[idx[q] for q in ((a-1,b),(a,b-1)) if q in idx] for (a,b) in P]
    f=[0]*len(P)
    def rec(t):
        if t==len(P):
            yield tuple(f[idx[(1,b)]] for b in range(1,(k[1] if len(k)>1 else 0)+1)); return
        lo=max([f[s] for s in preds[t]],default=0)+1
        for v in range(lo,M+1):
            if P[t][0]==P[t][1] and v%2: continue
            f[t]=v; yield from rec(t+1)
    yield from rec(0)
def H_at(k,n,x):
    d=sum(k); pts=list(range(0,d+1)); vals=[chiC(set(k),n,p)[0] for p in pts]
    s=Fr(0)
    for i,p in enumerate(pts):
        t=Fr(vals[i])
        for l,q in enumerate(pts):
            if l!=i: t*=Fr(x-q,p-q)
        s+=t
    return s
nmax=int(sys.argv[1]); A=int(sys.argv[2])
ok=bad=pos_ok=pos_bad=0
for L in range(1,nmax+1):
  for c in itertools.combinations(range(1,nmax+1),L):
    k=tuple(sorted(c,reverse=True)); n1=k[0]
    for a in range(1,A+1):
        M=2*a-1; tot=Fr(0)
        for u in rows_maps(k,M):
            p=Phi_def(n1,M,u)
            if M>=n1-1:
                if Phi_pos(n1,M,u)==p and p>0: pos_ok+=1
                else: pos_bad+=1; print('POSFAIL',k,a,u,p)
            else:
                if p==Fr((-1)**n1,2**n1): pos_ok+=1
                else: pos_bad+=1; print('SMALLFAIL',k,a,u,p)
            tot+=p
        h=(-1)**sum(k)*H_at(k,nmax,Fr(-2*a-1,2))
        if h==tot: ok+=1
        else: bad+=1; print('FAIL',k,a,h,tot)
print('n',nmax,'s_{a+1/2} = sum_g Phi vs Demazure:',ok,'ok',bad,'bad;  pointwise Phi lemma:',pos_ok,'ok',pos_bad,'bad')
