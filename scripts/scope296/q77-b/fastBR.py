"""dim B = sum_{D,D'} dim Hom_S(V x_D, V x_D'),  V x_D ~= V e_D (Norton x_D = pibar_{w0 D} pi_{w0 D^c}).
S acts through word maps phi_{mu,w}.  Linear algebra mod P (nmod_mat)."""
import sys, itertools, time, numpy as np
from collections import defaultdict
from flint import nmod_mat
from split import asc, move, act_pi, act_pibar, w0word
import os, resource
P=int(os.environ.get('P',str(2**31-1)))
MODE=os.environ.get('MODE','gen')
def isgen(mu,nu):
    d=[b-a for a,b in zip(mu,nu)]
    nz=[i for i,x in enumerate(d) if x]
    return len(nz)==2 and nz[1]==nz[0]+1 and abs(d[nz[0]])==1
NPROC=int(os.environ.get('NPROC','4'))
def mmod(A,B):
    B=B%P; lo=B&0xFFFF; hi=B>>16
    return (A@lo + ((A@hi)%P)*65536)%P
def mem(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss//1024
def nullity(rows,nc):
    if nc==0: return 0
    if not rows: return nc
    M=nmod_mat(len(rows),nc,[x%P for r in rows for x in r],P); return nc-M.rank()
def nullity2(blocks,nc):
    if nc==0: return 0
    if not blocks: return nc
    acc=None; buf=[]; cnt=0
    def flush(acc,buf):
        R=np.vstack(buf) if acc is None else np.vstack([acc]+buf)
        M=nmod_mat(R.shape[0],nc,[int(x) for x in R.ravel()],P); E,r=M.rref()
        return np.array(E.entries(),dtype=np.int64).reshape(R.shape)[:r]
    for blk in blocks:
        buf.append(blk); cnt+=blk.shape[0]
        if cnt>max(2000,2*nc):
            acc=flush(acc,buf); buf=[]; cnt=0
    if buf: acc=flush(acc,buf)
    return nc-acc.shape[0]
def run(n,N,want_perD=False):
    t0=time.time()
    W=list(itertools.product(range(1,n+1),repeat=N))
    wts=defaultdict(list)
    for w in W: wts[tuple(w.count(c) for c in range(1,n+1))].append(w)
    mus=sorted(wts)
    path={}
    for mu in mus:
        v0=tuple(c for c in range(1,n+1) for _ in range(mu[c-1]))
        Pp={v0:()}; fr=[v0]
        while fr:
            nf=[]
            for v in fr:
                for i in range(1,N):
                    u=move(v,i)
                    if u not in Pp: Pp[u]=Pp[v]+(i,); nf.append(u)
            fr=nf
        path[mu]=(v0,Pp)
    Smaps=[]  # (mu,nu,w)
    for mu in mus:
        S=asc(path[mu][0])
        for nu in mus:
            for w in wts[nu]:
                if asc(w)<=S and (MODE=='all' or isgen(mu,nu)): Smaps.append((mu,nu,w))
    print(n,N,'P',P,'MODE',MODE,'#constraint S-maps',len(Smaps),flush=True)
    def phi(mu,w,vec):
        Pp=path[mu][1]; out=defaultdict(int)
        for v,c in vec.items():
            u=w
            for i in Pp[v]: u=move(u,i)
            out[u]+=c
        return out
    Ds=sorted({asc(w) for w in W},key=lambda s:(len(s),sorted(s)))
    # basis of M_mu x_D: echelon in word coordinates
    def xD(vec,D):
        Dc=[j for j in range(1,N) if j not in D]
        for i in w0word(D): vec=act_pibar(vec,i)
        for i in w0word(Dc): vec=act_pi(vec,i)
        return vec
    Bas={}  # (D,mu) -> (rows as dict list, index of words, solver)
    for D in Ds:
        for mu in mus:
            ws=wts[mu]; ix={w:k for k,w in enumerate(ws)}
            vs=[xD({w:1},D) for w in ws]
            A=np.array([[v.get(w,0) for w in ws] for v in vs],dtype=np.int64)%P
            M=nmod_mat(A.shape[0],A.shape[1],[int(x) for x in A.ravel()],P)
            E,r=M.rref()
            ndw=sum(1 for w in ws if asc(w)==D)
            assert r==ndw, ('Norton dim mismatch',D,mu,r,ndw)
            if r==0: continue
            E=np.array(E.entries(),dtype=np.int64).reshape(A.shape)[:r]
            piv=[int(np.nonzero(e)[0][0]) for e in E]
            Bas[(D,mu)]=(E,ix,piv)
    print('bases done',time.time()-t0,flush=True)
    # S-action matrices on V x_D : for each smap and D, matrix (r_nu x r_mu) : coords of phi(basis) in basis of nu
    def coords(D,nu,vec):
        E,ix,piv=Bas[(D,nu)]
        v=np.zeros(len(ix),dtype=np.int64)
        for w,c in vec.items(): v[ix[w]]=c%P
        return v[piv]%P   # valid since vec in rowspace of rref E
    Act={D:{} for D in Ds}  # D -> {smap index: matrix r_nu x r_mu}
    for D in Ds:
        for si,(mu,nu,w) in enumerate(Smaps):
            if (D,mu) not in Bas or (D,nu) not in Bas: continue
            E,ix,piv=Bas[(D,mu)]; ws_mu=wts[mu]
            cols=[]
            for e in E:
                vec={ws_mu[k]:int(c) for k,c in enumerate(e) if c}
                cols.append(coords(D,nu,phi(mu,w,vec)))
            Mt=np.array(cols,dtype=np.int64).T
            if Mt.any(): Act[D][si]=Mt
    print('actions done',time.time()-t0,flush=True)
    rk=lambda D,mu: Bas[(D,mu)][0].shape[0] if (D,mu) in Bas else 0
    global _job
    def _job(D):
        rng=np.random.default_rng(hash((tuple(sorted(D)),P))%2**32)
        rowsum=0
        for D2 in Ds:
            sup=[mu for mu in mus if rk(D,mu) and rk(D2,mu)]
            if not sup: continue
            off={}; nv=0
            for mu in sup: off[mu]=nv; nv+=rk(D2,mu)*rk(D,mu)
            K=nv+8
            Mat=np.zeros((K,nv),dtype=np.int64)
            for si in set(Act[D])|set(Act[D2]):
                mu,nu,w=Smaps[si]
                M1=Act[D].get(si); M2=Act[D2].get(si)
                r1m,r2n=rk(D,mu),rk(D2,nu)
                if r1m==0 or r2n==0: continue
                U=rng.integers(0,P,size=(K,r2n),dtype=np.int64)
                V=rng.integers(0,P,size=(K,r1m),dtype=np.int64)
                if M1 is not None and nu in off:
                    r1n=rk(D,nu)
                    VM=mmod(V,M1.T)          # K x r1n
                    blk=(U[:,:,None]*VM[:,None,:])%P   # K x r2n x r1n
                    sl=slice(off[nu],off[nu]+r2n*r1n)
                    Mat[:,sl]=(Mat[:,sl]+blk.reshape(K,-1))%P
                if M2 is not None and mu in off:
                    r2m=rk(D2,mu)
                    UM=mmod(U,M2)            # K x r2m
                    blk=(UM[:,:,None]*V[:,None,:])%P   # K x r2m x r1m
                    sl=slice(off[mu],off[mu]+r2m*r1m)
                    Mat[:,sl]=(Mat[:,sl]-blk.reshape(K,-1))%P
            r=nmod_mat(K,nv,Mat.ravel().tolist(),P).rank()
            rowsum+=nv-r
        return D,rowsum
    import multiprocessing as mp
    tot=0; perD={}
    with mp.get_context('fork').Pool(NPROC) as pool:
        for D,rowsum in pool.imap_unordered(_job,sorted(Ds,key=lambda D:-sum(rk(D,m) for m in mus))):
            perD[D]=rowsum; tot+=rowsum
            print('  D',sorted(D),'row-sum',rowsum,'running',tot,'t',round(time.time()-t0,1),'MB',mem(),flush=True)
    print(n,N,'dim B_G =',tot,' time',round(time.time()-t0,1),'peakMB(main)',mem(),'children',resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss//1024,flush=True)
    if want_perD:
        for D in Ds: print('  D',sorted(D),'dim e_D B',perD[D])
    return tot,perD
if __name__=='__main__':
    for a in sys.argv[1:]:
        n,N=map(int,a.split(',')); run(n,N,want_perD=False)
