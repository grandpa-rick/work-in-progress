"""Sum certifier. A node = dict key->int, key=(T,lam,th,mu,thp,nu,psi,d,e) at level n with (a,b) and c-range.
Certified iff iterated leading coefficient (signed count at n=0) != 0."""
import sys, itertools
from functools import lru_cache
sys.path.insert(0,'../prove247')
from mp import subsets, geq, EA, EB, EC
U={'A':lambda n:2*n-4,'B':lambda n:2*n-2,'C':lambda n:2*n-2}
NV=None
LIMS=(0,1)
def supp(T,n,lam,mu,nu,a,b):
    if T=='A': return EA(n,a,b,lam,mu,nu)
    if T=='B': return EB(n,a,lam,mu,nu)
    return EC(n,lam,mu,nu)
def comp(n,nu): return tuple(k for k in range(1,n+1) if k not in nu)
def addW(e,upto,k=1):   # multiply by (prod_{j<upto} w_j)^k
    e=list(e)
    for j in range(1,upto): e[j]+=k
    return tuple(e)
def adde(e,j,k):
    e=list(e); e[j]+=k; return tuple(e)
def cexp(I):
    """c(I)=prod_{i in I}(-1)^i prod_{k<=i} w_k : (sign, exps)"""
    s=sum(I)%2; ex=[0]*(NV+1)
    for i in I:
        for k in range(1,i+1): ex[k]+=1
    return s,ex
def psi_opts(lam,th,n,lim):
    """options for Psi_{lam,th}(M), n in M. (order, sign, scalar e-delta (dict), ws, newlam, newth)"""
    r=len(lam)
    if lam and lam[-1]==n:
        if lim==0: return [(th[-1]-n-1,n%2,{k:1 for k in range(1,n)},-2,lam[:-1],th[:-1])]
        return [(th[-1],0,{},0,lam[:-1],th[:-1])]
    if lam==(1,) and th==(-1,) and lim==1: return [(-1,0,{},-2,(),())]
    if lam==(1,) and th==(1,) and lim==0: return [(-1,1,{},0,(),())]
    out=[]
    if lim==0:
        out.append((-r,r%2,{},0,lam,tuple(t-1 for t in th)))
        for j in range(r):
            J=j+1; o=th[j]-lam[j]-1-(r-J)
            out.append((o,(lam[j]+r-J)%2,{k:1 for k in range(1,lam[j])},-2,lam[:j]+lam[j+1:],th[:j]+tuple(t+1 for t in th[j+1:])))
    else:
        out.append((-r,r%2,{},0,lam,tuple(t+1 for t in th)))
        for j in range(r):
            J=j+1; o=th[j]-(r-J)
            out.append((o,(r-J)%2,{},0,lam[:j]+lam[j+1:],th[:j]+tuple(t+1 for t in th[j+1:])))
    return out
def apply_scalar(e,dd):
    e=list(e)
    for k,v in dd.items(): e[k]+=v
    return tuple(e)
def summand_opts(T,n,a,b,key,lim):
    """list of (const_order, slope, sign, newkey) ; newkey at level n-1 (e not yet cleaned of index n)"""
    T_,lam,th,mu,thp,nu,psi,d,e=key
    en=e[n]; res=[]
    if n not in nu:
        res.append((psi[-1]+en,0,0,(T,lam,th,mu,thp,nu,psi[:-1],d,e)))
        return res
    E=supp(T,n,lam,mu,nu,a,b)
    Ein=any(n in S for S in E); Eout=any(n not in S for S in E)
    if Ein:
        if lim==0: lo,ls,le,lw=1,n%2,{k:-1 for k in range(1,n)},0
        else: lo,ls,le,lw=-n,0,{},-2
        for A in psi_opts(lam,th,n,lim):
            for Bo in psi_opts(mu,thp,n,lim):
                o=lo+A[0]+Bo[0]; s=ls+A[1]+Bo[1]; ee=apply_scalar(apply_scalar(apply_scalar(e,le),A[2]),Bo[2])
                if T=='B':
                    o-=1; s+=1; ee=adde(ee,a,-1 if lim==0 else 1)
                nd=d+lw+A[3]+Bo[3]
                if T=='B' and lw: ee=adde(ee,a,lw)   # w_M^lw = w_a^lw w_S^lw
                res.append((o+en+d,1,s%2,(T,A[4],A[5],Bo[4],Bo[5],nu[:-1],psi,nd,ee)))
    if Eout:
        kc=comp(n,nu); km=nu[:-1]
        for D in psi_opts(kc,psi,n,lim):
            o,s,sc,ws,nk,npsi=D
            newnu=tuple(k for k in range(1,n) if k not in nk)
            s0,ex0=cexp(newnu); s1,ex1=cexp(km)
            ee=list(apply_scalar(e,sc))
            for k in range(1,NV+1): ee[k]+=ex0[k]-ex1[k]
            s+=s0+s1
            if lim==0:
                ee=addW(ee,n,-2); kk=ws+2; ee=addW(ee,n,kk); nd=d-kk; o=o+n
                if T=="B": ee=adde(ee,a,-kk)
            else:
                s+=n; ee=addW(ee,n,-1); ee=addW(ee,n,ws); nd=d-ws; o=o-1
                if T=="B": ee=adde(ee,a,-ws)
            res.append((o+en,0,s%2,(T,lam,th,mu,thp,newnu,npsi,nd,tuple(ee))))
    return res
def lower_side(labs,th,target):
    """rewrite Psi_{labs,th} in labels target (target<=labs positionwise): list of (sign, edelta dict, newth)"""
    out=[(0,{},tuple(th))]
    for i,(L,t) in enumerate(zip(labs,target)):
        if L==t: continue
        new=[]
        rng=list(range(t,L))
        for s,dd,th_ in out:
            for k in range(len(rng)+1):
                for Tset in itertools.combinations(rng,k):
                    d2=dict(dd)
                    for x in Tset: d2[x]=d2.get(x,0)+1
                    th2=list(th_); th2[i]-=k
                    new.append(((s+k)%2,d2,tuple(th2)))
        out=new
    return out
def canon(S,n):
    """S: dict key->coef at level n. lower labels per side to positionwise min among same-length."""
    for side in ('lam','mu','dual'):
        mins={}
        for key in S:
            T,lam,th,mu,thp,nu,psi,d,e=key
            labs = lam if side=='lam' else mu if side=='mu' else comp(n,nu)
            L=len(labs)
            mins[L]=labs if L not in mins else tuple(min(x,y) for x,y in zip(mins[L],labs))
        S2={}
        for key,cf in S.items():
            T,lam,th,mu,thp,nu,psi,d,e=key
            labs,tw = (lam,th) if side=='lam' else (mu,thp) if side=='mu' else (comp(n,nu),psi)
            tgt=mins[len(labs)]
            for s,dd,nth in lower_side(labs,tw,tgt):
                ee=apply_scalar(e,dd); sg=cf*(-1 if s else 1)
                if side=='lam': nk=(T,tgt,nth,mu,thp,nu,psi,d,ee)
                elif side=='mu': nk=(T,lam,th,tgt,nth,nu,psi,d,ee)
                else:
                    newnu=tuple(k for k in range(1,n+1) if k not in tgt)
                    s0,ex0=cexp(newnu); s1,ex1=cexp(nu)   # c(nu) stays; ell_{newnu} has c(newnu): multiply by c(newnu)/c(nu)
                    ee=list(ee)
                    for k in range(1,NV+1): ee[k]+=ex0[k]-ex1[k]
                    sg*= (-1)**((s0+s1)%2)
                    nk=(T,lam,th,mu,thp,newnu,nth,d,tuple(ee))
                S2[nk]=S2.get(nk,0)+sg
        S={k:v for k,v in S2.items() if v!=0}
    return S
def clean(S,n,a,b):
    return {k:v for k,v in S.items() if v!=0 and supp(k[0],n,k[1],k[3],k[5],a,b)}
STATS={'eqsplit':0}
def certify(S,n,a,b,cr,depth=0):
    S=clean(S,n,a,b)
    if not S: return False
    if n==0:
        return sum(S.values())!=0
    T=next(iter(S))[0]
    if T=='A' and n==b:
        S2={}
        for k,v in S.items():
            T_,lam,th,mu,thp,nu,psi,d,e=k
            if n not in nu: continue
            e=list(e); e[a]-=e[b]; e[b]=0
            nk=('B',lam,th,mu,thp,nu[:-1],psi,d,tuple(e)); S2[nk]=S2.get(nk,0)+v
        return certify(S2,n-1,a,None,cr,depth+1)
    if T=='B' and n==a:
        groups={}
        for k,v in S.items():
            T_,lam,th,mu,thp,nu,psi,d,e=k
            if n not in nu: continue
            ea=e[a]; e=list(e); e[a]=0
            nk=('C',lam,th,mu,thp,nu[:-1],psi,d-2,tuple(e))
            g=groups.setdefault(ea,{}); g[nk]=g.get(nk,0)+v
        for ea in sorted(groups):
            if certify(groups[ea],n-1,None,None,cr,depth+1): return True
        return False
    for lim in LIMS:
        opts=[]
        for k,v in S.items():
            for o,sl,sg,nk in summand_opts(T,n,a,b,k,lim):
                opts.append((o,sl,v*(-1 if sg else 1),nk))
        if not opts: continue
        A0=[o for o,sl,_,_ in opts if sl==0]; A1=[o for o,sl,_,_ in opts if sl==1]
        crs=[]
        if cr[0]=='eq':
            val=lambda o,sl:o+sl*cr[1]
            vals=[val(o,sl) for o,sl,_,_ in opts]; ext=min(vals) if lim==0 else max(vals)
            crs=[(cr,[x for x in opts if val(x[0],x[1])==ext])]
        else:
            cm=cr[1]
            if not A0 or not A1:
                vs=A0 or A1; ext=min(vs) if lim==0 else max(vs)
                crs=[(cr,[x for x in opts if x[0]==ext])]
            else:
                if lim==0:
                    m0,m1=min(A0),min(A1)
                    if m1+cm>m0: crs=[(cr,[x for x in opts if x[1]==0 and x[0]==m0])]
                    else:
                        top=m0-m1; STATS['eqsplit']+=1
                        crs=[(('eq',c),None) for c in range(cm,top+1)]+[(('ge',top+1),[x for x in opts if x[1]==0 and x[0]==m0])]
                else:
                    M0,M1=max(A0),max(A1)
                    if M1+cm>M0: crs=[(cr,[x for x in opts if x[1]==1 and x[0]==M1])]
                    else:
                        top=M0-M1; STATS['eqsplit']+=1
                        crs=[(('eq',c),None) for c in range(cm,top+1)]+[(('ge',top+1),[x for x in opts if x[1]==1 and x[0]==M1])]
        ok=True
        for cr2,sel in crs:
            if sel is None:
                if not certify(S,n,a,b,cr2,depth): ok=False;break
                continue
            ch={}
            for o,sl,cf,nk in sel:
                e=list(nk[8]); e[n]=0; nk=nk[:8]+(tuple(e),)
                ch[nk]=ch.get(nk,0)+cf
            ch={k:v for k,v in ch.items() if v}
            ch=canon(ch,n-1)
            if not certify(ch,n-1,a,b,cr2,depth+1): ok=False;break
        if ok: return True
    return False
def start(T,N,lam,mu,nu,a,b,cr=None):
    global NV
    NV=N
    th=tuple(lam[j]-2*(j+1) for j in range(len(lam))); thp=tuple(mu[j]-2*(j+1) for j in range(len(mu)))
    psi=tuple(0 for _ in comp(N,nu))
    d=2*len(lam)+2*len(mu)
    if cr is None: cr=('ge',U[T](N)+1)
    key=(T,lam,th,mu,thp,nu,psi,d,tuple([0]*(N+1)))
    return certify({key:1},N,a,b,cr)
def allA(N):
    for a,b in itertools.combinations(range(1,N+1),2):
        for nu in subsets(N):
            for lam in subsets(N):
                for mu in subsets(N):
                    if mu<lam: continue
                    if EA(N,a,b,lam,mu,nu): yield (lam,mu,nu,a,b)
