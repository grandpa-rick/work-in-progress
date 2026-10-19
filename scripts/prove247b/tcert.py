"""Twisted route certifier. Node: (T,n,lam,th,mu,thp,nu,psi,a,b,c) ; c = effective weight lower bound.
X = sum_M ell_{nu,psi}(e_M) w_M^c Psi_{lam,th}(M) Psi_{mu,thp}(M) [B: /g_a ; A: residue on H_ab]."""
import sys, itertools
from functools import lru_cache
sys.path.insert(0,'../prove247')
from mp import subsets, geq, EA, EB, EC
U={'A':lambda n:2*n-4,'B':lambda n:2*n-2,'C':lambda n:2*n-2}
FAIL=[]
def supp(T,n,lam,mu,nu,a,b):
    if T=='A': return EA(n,a,b,lam,mu,nu)
    if T=='B': return EB(n,a,lam,mu,nu)
    return EC(n,lam,mu,nu)
def opts0(lam,th):
    """u->0 options for Psi_{lam,th}(M), n in M, lam subset [n-1]: list of (order, weightshift, newlam, newth)"""
    if lam==(1,) and th==(1,): return [(-1,0,(),())]
    r=len(lam); out=[(-r,0,lam,tuple(t-1 for t in th))]
    for j in range(r):
        o=th[j]-lam[j]-1-(r-1-j)
        nl=lam[:j]+lam[j+1:]; nt=th[:j]+tuple(t+1 for t in th[j+1:])
        out.append((o,-2,nl,nt))
    return out
def optsI(lam,th):
    if lam==(1,) and th==(-1,): return [(-1,-2,(),())]
    r=len(lam); out=[(-r,0,lam,tuple(t+1 for t in th))]
    for j in range(r):
        o=th[j]-(r-1-j)
        nl=lam[:j]+lam[j+1:]; nt=th[:j]+tuple(t+1 for t in th[j+1:])
        out.append((o,0,nl,nt))
    return out
def uniq_min(op):
    m=min(o[0] for o in op); w=[o for o in op if o[0]==m]
    return w[0] if len(w)==1 else None
def uniq_max(op):
    m=max(o[0] for o in op); w=[o for o in op if o[0]==m]
    return w[0] if len(w)==1 else None
def comp(n,nu): return tuple(k for k in range(1,n+1) if k not in nu)
def fside(lam,th,n,lim):
    """leading of Psi_{lam,th}(M) for n in M. returns (ord,wshift,newlam,newth) or None if tie"""
    if lam and lam[-1]==n:   # forced
        if lim==0: return (th[-1]-n-1,-2,lam[:-1],th[:-1])
        return (th[-1],0,lam[:-1],th[:-1])
    return uniq_min(opts0(lam,th)) if lim==0 else uniq_max(optsI(lam,th))
def in_group(T,n,lam,th,mu,thp,nu,psi,c,lim):
    """leading of in-group (n in M), n in nu, nu^c subset [n-1]. returns (ord, child) or None"""
    A=fside(lam,th,n,lim); Bm=fside(mu,thp,n,lim)
    if A is None or Bm is None: return None
    if lim==0:
        o=1+c+A[0]+Bm[0]+(-1 if T=='B' else 0); cn=c+A[1]+Bm[1]
    else:
        o=-n+c+A[0]+Bm[0]+(-1 if T=='B' else 0); cn=c-2+A[1]+Bm[1]
    return o,(A[2],A[3],Bm[2],Bm[3],nu[:-1],psi,cn)
def out_group(T,n,lam,th,mu,thp,nu,psi,c,lim):
    kc=comp(n,nu)   # subset [n-1]
    op = opts0(kc,psi) if lim==0 else optsI(kc,psi)
    L = uniq_min(op) if lim==0 else uniq_max(op)
    if L is None: return None
    o,ws,nk,npsi=L
    if lim==0: o=o+n; cn=c-(ws+2)
    else: o=o-1; cn=c-ws
    newnu=tuple(k for k in range(1,n) if k not in nk)
    return o,(lam,th,mu,thp,newnu,npsi,cn)
@lru_cache(None)
def cert(T,n,lam,th,mu,thp,nu,psi,a,b,c):
    E=supp(T,n,lam,mu,nu,a,b)
    if not E: return False   # zero instance: cannot be a nonzero coefficient
    if T=='C' and n==0: return True
    if T=='A' and n==b: return cert('B',n-1,lam,th,mu,thp,nu[:-1],psi,a,None,c)
    if T=='B' and n==a: return cert('C',n-1,lam,th,mu,thp,nu[:-1],psi,None,None,c-2)
    if n not in nu:
        return cert(T,n-1,lam,th,mu,thp,nu,psi[:-1],a,b,c)
    Ein=[S for S in E if n in S]; Eout=[S for S in E if n not in S]
    tried=False
    for lim in (0,1):
        I=in_group(T,n,lam,th,mu,thp,nu,psi,c,lim) if Ein else None
        O=out_group(T,n,lam,th,mu,thp,nu,psi,c,lim) if Eout else None
        if Ein and Eout:
            if I is None or O is None:
                # dominant group might still be clean; need the other group's termwise bound
                pass
            if lim==0:
                # out dominates iff out order < in order ; in order is termwise bound min (use even if tie)
                if O is not None:
                    io=in_bound(T,n,lam,th,mu,thp,c,0)
                    if O[0]<io:
                        tried=True
                        if cert(T,n-1,*O[1][:6],a,b,O[1][6]): return True
                if I is not None:
                    ob=out_bound(n,nu,psi,0)
                    if I[0]<ob:
                        tried=True
                        if cert(T,n-1,*I[1][:6],a,b,I[1][6]): return True
            else:
                if I is not None:
                    ob=out_bound(n,nu,psi,1)
                    if I[0]>ob:
                        tried=True
                        if cert(T,n-1,*I[1][:6],a,b,I[1][6]): return True
                if O is not None:
                    io=in_bound(T,n,lam,th,mu,thp,c,1)
                    if O[0]>io:
                        tried=True
                        if cert(T,n-1,*O[1][:6],a,b,O[1][6]): return True
        else:
            G = I if Ein else O
            if G is not None:
                tried=True
                if cert(T,n-1,*G[1][:6],a,b,G[1][6]): return True
    if not tried: FAIL.append((T,n,lam,th,mu,thp,nu,psi,a,b,c))
    return False
def fbound(lam,th,n,lim):
    if lam and lam[-1]==n: return fside(lam,th,n,lim)[0]
    op=opts0(lam,th) if lim==0 else optsI(lam,th)
    return min(o[0] for o in op) if lim==0 else max(o[0] for o in op)
def in_bound(T,n,lam,th,mu,thp,c,lim):
    s=fbound(lam,th,n,lim)+fbound(mu,thp,n,lim)+(-1 if T=='B' else 0)
    return 1+c+s if lim==0 else -n+c+s
def out_bound(n,nu,psi,lim):
    kc=comp(n,nu); op=opts0(kc,psi) if lim==0 else optsI(kc,psi)
    return (min(o[0] for o in op)+n) if lim==0 else (max(o[0] for o in op)-1)
def start(T,N,lam,mu,nu,a,b):
    th=tuple(lam[j]-2*(j+1) for j in range(len(lam))); thp=tuple(mu[j]-2*(j+1) for j in range(len(mu)))
    psi=tuple(0 for _ in comp(N,nu)); c=U[T](N)+1+2*len(lam)+2*len(mu)
    return cert(T,N,lam,th,mu,thp,nu,psi,a,b,c)
def allA(N):
    for a,b in itertools.combinations(range(1,N+1),2):
        for nu in subsets(N):
            for lam in subsets(N):
                for mu in subsets(N):
                    if mu<lam: continue
                    if EA(N,a,b,lam,mu,nu): yield (lam,mu,nu,a,b)
if __name__=='__main__':
  import collections
  for N in range(3,int(sys.argv[1])+1):
    FAIL.clear(); cert.cache_clear(); tot=bad=0; badl=[]
    for lam,mu,nu,a,b in allA(N):
        tot+=1
        if not start('A',N,lam,mu,nu,a,b): bad+=1; badl.append((lam,mu,nu,a,b))
    print(N,'allowed',tot,'uncert',bad,'nofail-route nodes',len(set(FAIL)),collections.Counter((f[0],f[1]) for f in set(FAIL)))
    print('  ',badl[:5]); sys.stdout.flush()
