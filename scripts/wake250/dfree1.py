# Task 1+3: compare Willems classes with Psi_{nu,0}(M; w or 1/w) under candidate dictionaries.
import sys, random
from fractions import Fraction as Fr
sys.path.insert(0,'../prove249c')
from eng import Psi, subsets, leq
from kschub import Cn
n=int(sys.argv[1]); seed=int(sys.argv[2]) if len(sys.argv)>2 else 0
if seed==0: w={i:Fr([5,7,11,13,17][i-1]) for i in range(1,n+1)}
else:
    random.seed(seed); w={i:Fr(random.choice([5,7,11,13,17,19,23,29,31]),random.choice([37,41,43,47,53,59])) for i in range(1,n+1)}
C=Cn(n); S=[tuple(x) for x in subsets(n)]
def sv(M,rev=True):  # sign vector: -1 at positions n+1-k (rev) or k
    v=[1]*n
    for k in M: v[(n-k) if rev else (k-1)]=-1
    return tuple(v)
comp=lambda M: tuple(k for k in range(1,n+1) if k not in M)
print('n=%d seed=%d'%(n,seed))
print('length check l(w_nu)=|nu| (rev positions):',all(C.len[C.rep[sv(v)]]==sum(v) for v in S))
print('Bruhat check: support of Willems class = {M: leq(nu,M)}? (computed below)')
inv=lambda d:{k:1/d[k] for k in d}
dicts={}
for rev in (True,False):
    for tb in ('w','1/w'):
        def tmap(rev=rev,tb=tb):
            base=w if tb=='w' else inv(w)
            return [base[n-p] if rev else base[p+1] for p in range(n)]  # p=0..n-1 -> position p+1 ; rev: w_{n+1-(p+1)}=w_{n-p}
        dicts[(rev,tb)]=tmap()
Pt={(b,v,M):Psi(list(v),[0]*len(v),list(M),w if b=='w' else inv(w)) for b in ('w','1/w') for v in S for M in S}
res=[]
for (rev,tb),t in dicts.items():
    Wl={(v,M):C.willems(C.rep[sv(v,rev)],C.rep[sv(M,rev)],t) for v in S for M in S}
    supp_ok=all((Wl[v,M]!=0)==leq(list(v),list(M)) for v in S for M in S)
    for b in ('w','1/w'):
        exact=all(Wl[v,M]==Pt[b,v,M] for v in S for M in S)
        prop=0
        for v in S:
            r=set()
            bad=False
            for M in S:
                a,c=Wl[v,M],Pt[b,v,M]
                if a==0 and c==0: continue
                if a==0 or c==0: bad=True;break
                r.add(a/c)
            if not bad and len(r)==1: prop+=1
        line='positions=%s t=%-3s | Willems support==dominance:%s | Psi base=%-3s exact=%s  prop-per-nu %d/%d'%('rev' if rev else 'id',tb,supp_ok,b,exact,prop,len(S))
        print(line); res.append(((rev,tb,b),exact,prop))
