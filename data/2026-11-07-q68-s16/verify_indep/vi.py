"""Independent verifier (Wake 282). Fresh code; no imports from prove280/prove281/wake282.
Conventions: one-line word w=[w(1..n)]; P = RSK insertion tableau of the word, Q = recording.
Left mult by s_i swaps VALUES i,i+1. L(w) = {i : i+1 appears before i in the word} = Des P(w)
(Des T = {i : i+1 in a strictly lower row than i}). Left cell = {Q fixed}."""
import sys, random
from itertools import combinations

def syt_rowwords(lam):
    """yield row-word r[0..n-1]: r[k] = row of value k+1."""
    n=sum(lam); L=len(lam); fill=[0]*L; r=[0]*n
    def rec(k):
        if k==n: yield tuple(r); return
        for i in range(L):
            if fill[i]<lam[i] and (i==0 or fill[i-1]>fill[i]):
                fill[i]+=1; r[k]=i
                yield from rec(k+1)
                fill[i]-=1
    yield from rec(0)

def rw_to_tab(rw,lam): return [[k+1 for k in range(len(rw)) if rw[k]==i] for i in range(len(lam))]
def tab_des(T):
    row={v:i for i,R in enumerate(T) for v in R}; n=len(row)
    return frozenset(i for i in range(1,n) if row[i+1]>row[i])
def rw_desmask(rw):
    m=0
    for i in range(1,len(rw)):
        if rw[i]>rw[i-1]: m|=1<<i   # bit i <-> descent i
    return m
def interior(A,n):
    return frozenset(i for i in A if all((j in A) for j in (i-1,i+1) if 1<=j<=n-1))

def rsk(w):
    P=[];Q=[]
    for pos,a in enumerate(w,1):
        r=0
        while True:
            if r==len(P): P.append([a]); Q.append([pos]); break
            row=P[r]; j=next((k for k,b in enumerate(row) if b>a),None)
            if j is None: row.append(a); Q[r].append(pos); break
            row[j],a=a,row[j]; r+=1
    return P,Q

def inv_rsk(P,Q):
    P=[list(r) for r in P]; Q=[list(r) for r in Q]; n=sum(map(len,P)); w=[0]*n
    for pos in range(n,0,-1):
        r=next(i for i,R in enumerate(Q) if R and R[-1]==pos)
        Q[r].pop(); a=P[r].pop()
        for rr in range(r-1,-1,-1):
            row=P[rr]; j=max(k for k,b in enumerate(row) if b<a)
            row[j],a=a,row[j]
        P=[R for R in P if R] ; Q=[R for R in Q if R]
        w[pos-1]=a
    return w

def length(w): return sum(1 for i in range(len(w)) for j in range(i+1,len(w)) if w[i]>w[j])
def Ldes(w):
    pos={v:i for i,v in enumerate(w)}; return frozenset(i for i in range(1,len(w)) if pos[i+1]<pos[i])
def lmul(i,w): return [i+1 if v==i else i if v==i+1 else v for v in w]

def bruhat_le(u,w):
    """rank-matrix criterion: u<=w iff for all i,j: #{a<=i: u(a)>=j} <= #{a<=i: w(a)>=j}."""
    n=len(u)
    cu=[0]*(n+2); cw=[0]*(n+2)
    for i in range(n):
        for j in range(1,u[i]+1): cu[j]+=1
        for j in range(1,w[i]+1): cw[j]+=1
        for j in range(1,n+1):
            if cu[j]>cw[j]: return False
    return True

def mu_cheap(a,b):
    """returns 0 if forced zero by parity/incomparability, 1 if gap-1 comparable, None otherwise"""
    la,lb=length(a),length(b)
    if la==lb: return 0
    lo,hi=(a,b) if la<lb else (b,a)
    if not bruhat_le(lo,hi): return 0
    if (abs(la-lb))%2==0: return 0
    if abs(la-lb)==1: return 1
    # F1: s in L(hi)\L(lo) and lo != s hi  => mu=0
    for s in Ldes(hi)-Ldes(lo):
        if lmul(s,hi)!=lo: return 0
    return None

def lstar(j,w):
    """left star wrt {s_j,s_{j+1}}; None if w not in D_L."""
    L=Ldes(w); t=(j in L, j+1 in L)
    if t[0]==t[1]: return None
    for k in (j,j+1):
        v=lmul(k,w); Lv=Ldes(v); tv=(j in Lv, j+1 in Lv)
        if tv[0]!=tv[1] and tv!=t: return v
    raise RuntimeError

def star_search(a,b,maxnodes=200000):
    """BFS over simultaneous left stars; mu~ invariant (Lemma S). return (value, path)."""
    from collections import deque
    n=len(a); start=(tuple(a),tuple(b)); seen={start:None}; dq=deque([start])
    while dq:
        p=dq.popleft(); m=mu_cheap(list(p[0]),list(p[1]))
        if m is not None:
            path=[];q=p
            while seen[q] is not None: q,j=seen[q]; path.append(j)
            return m,path[::-1]
        if len(seen)>maxnodes: break
        for j in range(1,n-1):
            u=lstar(j,list(p[0])); v=lstar(j,list(p[1]))
            if u is None or v is None: continue
            q=(tuple(u),tuple(v))
            if q not in seen: seen[q]=(p,j); dq.append(q)
    return None,None
