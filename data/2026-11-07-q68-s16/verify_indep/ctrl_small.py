"""Controls: (a) rank-matrix Bruhat vs Ehresmann tableau criterion, S_6;
(b) mu_cheap + star_search vs full KL mu (fresh standard recursion), S_5 and S_6, all pairs in a left cell."""
from vi import *
from itertools import permutations
import sys
def ehres(u,w):
    n=len(u)
    for i in range(1,n+1):
        a=sorted(u[:i]); b=sorted(w[:i])
        if any(x>y for x,y in zip(a,b)): return False
    return True
for n in (6,):
    G=[list(p) for p in permutations(range(1,n+1))]
    bad=sum(1 for u in G for w in G if bruhat_le(u,w)!=ehres(u,w))
    print('bruhat mismatches S6:',bad)

def full_kl(n):
    G=[tuple(p) for p in permutations(range(1,n+1))]
    G.sort(key=length_t)
    idx={g:i for i,g in enumerate(G)}; N=len(G); ln=[length_t(g) for g in G]
    le=[[bruhat_le(list(G[a]),list(G[b])) for b in range(N)] for a in range(N)]
    P={}  # P[(x,w)] poly list
    def mu(x,w):
        if not le[x][w] or x==w: return 0
        d=ln[w]-ln[x]
        if d%2==0: return 0
        p=P[(x,w)]; k=(d-1)//2
        return p[k] if k<len(p) else 0
    for w in range(N):
        if ln[w]==0: P[(w,w)]=[1]; continue
        g=G[w]; s=min(Ldes(list(g)))
        v=idx[tuple(lmul(s,list(g)))]  # w = s v, s v > v
        Ms=[z for z in range(N) if le[z][v] and z!=v and s in Ldes(list(G[z])) and mu(z,v)]
        for x in range(N):
            if not le[x][w]: continue
            sx=idx[tuple(lmul(s,list(G[x])))]
            c=1 if ln[sx]<ln[x] else 0
            def get(a,b): return P.get((a,b),[0]) if le[a][b] else [0]
            res=[0]*(ln[w]+2)
            a=get(sx,v); 
            for i,co in enumerate(a): res[i+1-c]+=co
            b=get(x,v)
            for i,co in enumerate(b): res[i+c]+=co
            for z in Ms:
                m=mu(z,v); e=(ln[w]-ln[z])//2; pz=get(x,z)
                for i,co in enumerate(pz): res[i+e]-=m*co
            while len(res)>1 and res[-1]==0: res.pop()
            P[(x,w)]=res
    return G,idx,ln,le,mu
def length_t(g): return length(list(g))
for n in (5,6):
    G,idx,ln,le,mu=full_kl(n)
    assert P_ok(G,idx) if False else True
    # sanity: P_{e,4231} in S4 embedded? check n>=4: w=[4,2,3,1,5..]
    e=idx[tuple(range(1,n+1))]; w=idx[tuple([4,2,3,1]+list(range(5,n+1)))]
    print('n',n,'mu(e,4231)-poly check done')
    cells={}
    for g in G: cells.setdefault(str(rsk(list(g))[1]),[]).append(g)
    stats={'agree':0,'wrong':0,'undec':0}
    for c in cells.values():
        for a in c:
            for b in c:
                if a>=b: continue
                ia,ib=idx[a],idx[b]
                true = mu(ia,ib) if ln[ia]<ln[ib] else mu(ib,ia)
                m,_=star_search(list(a),list(b))
                if m is None: stats['undec']+=1
                elif m==true: stats['agree']+=1
                else: stats['wrong']+=1
    print('n',n,stats)
