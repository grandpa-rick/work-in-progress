import itertools, sys, numpy as np, math, random
from mpq import S, subsets, geq, Fv, EA
from fractions import Fraction as Fr
def build(N,a,b,al,dl):
    w=[S.const(1+al[i]*dl) for i in range(N)]
    w[b-1]=S(0,[1/(1+al[a-1]*dl)]*6)
    subs=subsets(N); idx={M:i for i,M in enumerate(subs)}
    Fm={}
    for I in subs:
        for M in subs:
            if geq(M,I): Fm[(I,M)]=Fv(I,M,w)
    # inverse of F (rows I, cols M), F[I][M] nonzero only if M>=I : ell = f F^{-1}
    # solve per nu: ell_nu(e_M) via triangular recursion: coord of delta_M
    G={}
    order=sorted(subs,key=lambda T:(len(T),sum(T)))
    for M in subs:
        c={}
        for X in order:
            if not geq(X,M): continue
            v=S.const(1 if X==M else 0)
            for I in c:
                if geq(X,I): v=v-c[I]*Fm[(I,X)]
            c[X]=v/Fm[(X,X)]
        for X in c: G[(M,X)]=c[X]       # ell_X(e_M)
    return w,subs,Fm,G
def RA_all(N,a,b,lam,mu,nu,cs,sysd):
    w,subs,Fm,G=sysd; out=[]
    terms=[]
    for M in subs:
        if a in M and b in M and geq(nu,M) and (M,nu) in G:
            Mpp=tuple(m for m in M if m not in (a,b))
            g=G[(M,nu)]
            g=g*(Fm[(lam,Mpp)] if (lam,Mpp) in Fm else S.const(0))*(Fm[(mu,Mpp)] if (mu,Mpp) in Fm else S.const(0))
            terms.append((Mpp,g))
    for c in cs:
        tot=S.const(0)
        for Mpp,g in terms:
            t=g
            for m in Mpp: t=t*(w[m-1]**c)
            tot=tot+t
        out.append(tot.coef(-1))
    return out
if __name__=="__main__":
    N=int(sys.argv[1]); random.seed(int(sys.argv[2]))
    al=[Fr(random.randint(1,20),random.randint(1,20))*random.choice([1,-1]) for _ in range(N)]
    Uall=list(range(1,N+1)); subs=subsets(N)
    dim=N*(N+1)//2; Umax=2*N-4
    cs=list(range(Umax+1-2*(dim+2),Umax+2+2*(dim+2)))
    dls=[Fr(k,10**14) for k in (1,2,4)]
    pairs=eval(sys.argv[3]) if len(sys.argv)>3 else [(a,b) for a in Uall for b in Uall if a<b]
    normsys={dl:build(N,N-1,N,al,dl) for dl in dls}
    normv={dl:RA_all(N,N-1,N,(),(),tuple(Uall),[2*N-2],normsys[dl])[0] for dl in dls}
    for (a,b) in pairs:
     sy={dl:build(N,a,b,al,dl) for dl in dls}
     for nu in subs:
      for lam in subs:
       for mu in subs:
        if lam>mu: continue
        if not EA(N,a,b,lam,mu,nu): continue
        V={dl:[v/normv[dl] for v in RA_all(N,a,b,lam,mu,nu,cs+[Umax+7],sy[dl])] for dl in dls}
        v0=V[dls[0]][-1]; v1=V[dls[1]][-1]
        if v0==0: print('ZERO at U+7',(a,b),lam,mu,nu); continue
        d=round(math.log2(abs(float(v1/v0))))
        vals=[]
        for i in range(len(cs)):
            g=[V[dl][i]/dl**d for dl in dls]
            vals.append(float((8*g[0]-6*g[1]+g[2])/3))
        co=np.polyfit(cs,vals,dim+2); sc=max(abs(v) for v in vals)
        co=np.trim_zeros(np.where(abs(co)*max(abs(np.array(cs)))**np.arange(len(co))[::-1]<1e-7*sc,0,co),'f')
        r=np.roots(co) if len(co)>1 else []
        mr=max([x.real for x in r],default=None)
        realr=sorted(round(x.real,3) for x in r if abs(x.imag)<1e-6)
        intabove=[c for c in range(Umax+1,Umax+40) if abs(np.polyval(co,c))<1e-6*sc]
        print((a,b),lam,mu,nu,'d',d,'deg',len(co)-1,'maxRe',None if mr is None else round(mr,3),'real',realr,'intzeros>U',intabove,flush=True)
