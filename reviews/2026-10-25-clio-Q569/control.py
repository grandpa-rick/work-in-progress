"""Independent control (Rick 2026-10-25): hex-puzzles vs K-homology coproduct coefficients, n=5,k=2,
and hexR-puzzles vs d~ = d on transposes. Built only from Clio's printed PDF defs (sec 2, Conv 2, sec 3)
+ set-valued-tableau G_lambda. No code of hers used."""
import itertools, sys, numpy as np
from collections import defaultdict
from tiles import TILES, labelled, normalise, tri_pts
# ---------- Grothendieck polynomials via set-valued tableaux ----------
def svt_poly(lam, nvars):
    cells=[(i,j) for i,r in enumerate(lam) for j in range(r)]
    subsets=[frozenset(c) for k in range(1,nvars+1) for c in itertools.combinations(range(nvars),k)]
    P=defaultdict(int); T={}
    def rec(idx):
        if idx==len(cells):
            mon=[0]*nvars; sz=0
            for S in T.values():
                for x in S: mon[x]+=1
                sz+=len(S)
            P[tuple(mon)]+=(-1)**(sz-len(cells)); return
        i,j=cells[idx]
        for S in subsets:
            if j>0 and min(S)<max(T[(i,j-1)]): continue
            if i>0 and min(S)<=max(T[(i-1,j)]): continue
            T[(i,j)]=S; rec(idx+1); del T[(i,j)]
    rec(0); return {m:c for m,c in P.items() if c}
def parts_in(rows,cols):
    out=[]
    def rec(pre,mx,left):
        if left==0: out.append(tuple(p for p in pre if p)); return
        for v in range(mx,-1,-1): rec(pre+[v],v,left-1)
    rec([],cols,rows); return sorted(set(out),key=lambda p:(sum(p),p))
def conj(p): return tuple(sum(1 for x in p if x>j) for j in range(p[0])) if p else ()
_cache={}
def G(lam,N):
    if (lam,N) not in _cache: _cache[(lam,N)]=svt_poly(lam,N)
    return _cache[(lam,N)]
def coproduct(nu,N):
    """return dict (lam,mu)->coef with G_nu(x;y)=sum coef G_lam(x)G_mu(y), N vars each"""
    full=svt_poly(nu,2*N)
    if len(nu)>N: raise ValueError
    basis=[(l,m) for l in parts_in(N,nu[0] if nu else 0) for m in parts_in(N,nu[0] if nu else 0)]
    cols=[]; mons={}
    for l,m in basis:
        gl,gm=G(l,N),G(m,N); prod={}
        for a,ca in gl.items():
            for b,cb in gm.items(): prod[a+b]=prod.get(a+b,0)+ca*cb
        cols.append(prod)
        for k in prod: mons.setdefault(k,len(mons))
    for k in full: mons.setdefault(k,len(mons))
    A=np.zeros((len(mons),len(basis))); y=np.zeros(len(mons))
    for j,c in enumerate(cols):
        for k,v in c.items(): A[mons[k],j]=v
    for k,v in full.items(): y[mons[k]]=v
    sol=np.rint(np.linalg.lstsq(A,y,rcond=None)[0]).astype(int)
    assert np.array_equal(A@sol,y), 'expansion not exact'
    return {basis[j]:int(sol[j]) for j in range(len(basis)) if sol[j]}
# ---------- puzzles ----------
D=[(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
def sigma(lam,n,k):
    lam=list(lam)+[0]*(k-len(lam)); s=[0]*(n-k-lam[0])
    for i in range(k):
        s.append(1); s+= [0]*(lam[i]-(lam[i+1] if i+1<k else 0))
    assert len(s)==n and sum(s)==k, (lam,s); return s
def count_puzzles(tileset,n,lam,mu,nu,k):
    region={('U',a,b) for a in range(n) for b in range(n) if a+b<=n-1}|{('D',a,b) for a in range(n) for b in range(n) if a+b<=n-2}
    E={}
    sl,sm,sn=sigma(lam,n,k),sigma(mu,n,k),sigma(nu,n,k)
    for i in range(n):
        E[frozenset(((i,0),(i+1,0)))]=sn[i]
        E[frozenset(((0,i),(0,i+1)))]=sl[i]
        E[frozenset(((i,n-i),(i+1,n-i-1)))]=sm[i]
    order=sorted(region,key=lambda t:(t[2],t[1],t[0]=='D'))
    places=defaultdict(list)  # triangle -> list of (tris,edges)
    for tris,ed in tileset:
        for t0 in tris:
            for t in region:
                if t[0]!=t0[0]: continue
                da,db=t[1]-t0[1],t[2]-t0[2]
                tt=frozenset((x[0],x[1]+da,x[2]+db) for x in tris)
                if not tt<=region: continue
                ee=[(frozenset((p[0]+da,p[1]+db) for p in e),l) for e,l in ed]
                places[t].append((tt,ee))
    covered=set(); cnt=0
    def rec(i):
        nonlocal cnt
        while i<len(order) and order[i] in covered: i+=1
        if i==len(order): cnt+=1; return
        for tt,ee in places[order[i]]:
            if tt&covered: continue
            if any(E.get(e,l)!=l for e,l in ee): continue
            new=[e for e,l in ee if e not in E]
            for e,l in ee: E.setdefault(e,l)
            covered.update(tt); rec(i+1); covered.difference_update(tt)
            for e in new: del E[e]
    rec(0); return cnt
NT={nm:normalise(labelled(s)) for nm,s in TILES.items()}
KTW=[NT[x] for x in ['tri0U','tri0D','tri1U','tri1D','rh0','rh1','rh2']]
def sgn(e): return 1 if e%2==0 else -1
if __name__=='__main__':
    n,k=5,2
    # sanity: KTW = classical LR on degree-equal triples (k x (n-k))
    P=parts_in(k,n-k); bad=0; tot=0
    copro={nu:coproduct(nu,k) for nu in P}
    for l,m,nu in itertools.product(P,P,P):
        if sum(nu)!=sum(l)+sum(m): continue
        tot+=1; c=copro[nu].get((l,m),0)
        if count_puzzles(KTW,n,l,m,nu,k)!=c: bad+=1
    print('KTW vs LR (coproduct, top degree), n=5,k=2:',tot-bad,'/',tot)
    # hex control
    Ph=parts_in(k,n-k-1); ok=tot=nz=0; fails=[]
    for l,m,nu in itertools.product(Ph,Ph,Ph):
        if sum(nu)>sum(l)+sum(m): continue
        tot+=1
        d=sgn(sum(l)+sum(m)-sum(nu))*copro[nu].get((l,m),0)
        p=count_puzzles(KTW+[NT['hex']],n,l,m,nu,k)
        nz+= d!=0
        if p==d: ok+=1
        else: fails.append((l,m,nu,p,d))
    print('hex-puzzles vs d, n=5,k=2:',ok,'/',tot,'nonzero d:',nz,'fails',fails[:5])
    # hexR control: (k-1)x(n-k), d~ = d on transposes
    PR=parts_in(k-1,n-k); ok=tot=nz=0; fails=[]
    cop3={}
    for l,m,nu in itertools.product(PR,PR,PR):
        if sum(nu)>sum(l)+sum(m): continue
        tot+=1; nc=conj(nu)
        if nc not in cop3: cop3[nc]=coproduct(nc,max(1,len(nc)))
        d=sgn(sum(l)+sum(m)-sum(nu))*cop3[nc].get((conj(l),conj(m)),0)
        p=count_puzzles(KTW+[NT['hexR']],n,l,m,nu,k)
        nz+= d!=0
        if p==d: ok+=1
        else: fails.append((l,m,nu,p,d))
    print('hexR-puzzles vs d~, n=5,k=2:',ok,'/',tot,'nonzero:',nz,'fails',fails[:5])
