"""dim A = sum_D dim(gH), g = (b_1..b_r) in V^r, b_i basis of V x_D; gH by Krylov closure under pi_i (mod P)."""
import sys, itertools, time, numpy as np
from flint import nmod_mat
from split import asc, move, act_pi, act_pibar, w0word
import os
P=int(os.environ.get('P',str(2**31-1)))
class Ech:
    def __init__(s): s.rows=[]; s.piv=[]
    def add(s,v):
        v=v%P
        for r,p in zip(s.rows,s.piv):
            if v[p]: v=(v-v[p]*r)%P
        nz=np.nonzero(v)[0]
        if len(nz)==0: return False
        p=nz[0]; inv=pow(int(v[p]),P-2,P); v=(v*inv)%P
        s.rows.append(v); s.piv.append(p); return True
def run(n,N):
    t0=time.time()
    W=list(itertools.product(range(1,n+1),repeat=N)); widx={w:k for k,w in enumerate(W)}; d=len(W)
    gens=[np.array([widx[move(w,i)] for w in W]) for i in range(1,N)]
    Ds=sorted({asc(w) for w in W},key=lambda s:(len(s),sorted(s)))
    def xD(vec,D):
        Dc=[j for j in range(1,N) if j not in D]
        for i in w0word(D): vec=act_pibar(vec,i)
        for i in w0word(Dc): vec=act_pi(vec,i)
        return vec
    tot=0; PER={}
    for D in Ds:
        vs=[xD({w:1},D) for w in W]; vs=[v for v in vs if v]
        A=np.zeros((len(vs),d),dtype=np.int64)
        for r,v in enumerate(vs):
            for w,c in v.items(): A[r,widx[w]]=c%P
        E,r=nmod_mat(A.shape[0],d,[int(x) for x in A.ravel()],P).rref()
        E=np.array(E.entries(),dtype=np.int64).reshape(A.shape)[:r]
        # g as r x d matrix; action of pi_i: (b.pi_i)[f(w)] += b[w]
        def act(G,i):
            out=np.zeros_like(G); f=gens[i]
            for k in range(G.shape[0]): np.add.at(out[k],f,G[k])
            return out%P
        ech=Ech(); g=E.copy(); queue=[g]; ech.add(g.ravel())
        while queue:
            X=queue.pop()
            for i in range(N-1):
                Y=act(X,i)
                if ech.add(Y.ravel()): queue.append(Y)
        PER[tuple(sorted(D))]=len(ech.rows); tot+=len(ech.rows)
    print(n,N,'dim A =',tot,'time',round(time.time()-t0,1),flush=True)
    return tot,PER
if __name__=='__main__':
    for a in sys.argv[1:]:
        n,N=map(int,a.split(',')); run(n,N)
