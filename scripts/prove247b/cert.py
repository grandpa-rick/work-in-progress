import sys, itertools
from functools import lru_cache
sys.path.insert(0,'../prove247')
from mp import subsets, geq, EA, EB, EC
U={'A':lambda n:2*n-4,'B':lambda n:2*n-2,'C':lambda n:2*n-2}
def top(l): return l[-1] if l else None
def minus(l): return l[:-1]
def tclean(l): return len(l)==0 or l==(1,) or (l[-1]-1 not in l)
def nutil(nu,n):
    comp=[k for k in range(1,n+1) if k not in nu]; m=max(comp)
    return tuple(sorted([x for x in nu if x!=n]+[m])), m, comp
FAIL=[]
def supp(T,n,lam,mu,nu,a,b):
    if T=='A': return EA(n,a,b,lam,mu,nu)
    if T=='B': return EB(n,a,lam,mu,nu)
    return EC(n,lam,mu,nu)
@lru_cache(None)
def cert(T,n,lam,mu,nu,a,b,routes):
    """routes: set of allowed mixed routes ('0','inf'). returns True if certified nonzero for all c>U^T_n, given support nonempty."""
    E=supp(T,n,lam,mu,nu,a,b)
    if not E: return True   # nothing to prove
    if T=='C' and n==0: return True
    if T=='A' and n==b: return cert('B',n-1,lam,mu,minus(nu),a,None,routes)
    if T=='B' and n==a: return cert('C',n-1,lam,mu,minus(nu),None,None,routes)
    if n not in nu: return cert(T,n-1,lam,mu,nu,a,b,routes)
    if n in lam or n in mu: return cert(T,n-1,minus(lam),minus(mu),minus(nu),a,b,routes)
    Ein=[S for S in E if n in S]; Eout=[S for S in E if n not in S]
    delta=(lam!=())+(mu!=())
    if not Ein:  # stuck: shift (c-free)
        return cert(T,n-1,lam,mu,nutil(nu,n)[0],a,b,routes)
    if not Eout:
        return cert(T,n-1,minus(lam),minus(mu),minus(nu),a,b,routes)
    tried=False
    cmin=U[T](n)+1
    inord = cmin-delta if T=='B' else cmin+1-delta
    if '0' in routes:
        if len(nu)==n: o=n; ok=True; nn=tuple(range(1,n))
        else:
            nn,m,comp=nutil(nu,n); o=n-m-1
            ok = len(comp)==1 or comp[-2]<m-1
        if ok and o<inord:
            if cert(T,n-1,lam,mu,nn,a,b,routes): return True
            tried=True
    if 'inf' in routes and tclean(lam) and tclean(mu):
        indeg = cmin-n-(1 if T=='B' else 0)
        if indeg> -1:
            if cert(T,n-1,minus(lam),minus(mu),minus(nu),a,b,routes): return True
            tried=True
    if not tried: FAIL.append((T,n,lam,mu,nu,a,b))
    return False
def allA(N):
    for a,b in itertools.combinations(range(1,N+1),2):
        for nu in subsets(N):
            for lam in subsets(N):
                for mu in subsets(N):
                    if mu<lam: continue
                    if EA(N,a,b,lam,mu,nu): yield (lam,mu,nu,a,b)
if __name__=='__main__':
  for N in range(3,int(sys.argv[1])+1):
    for routes in (('inf',),('0',),('0','inf')):
        FAIL.clear(); cert.cache_clear()
        tot=bad=0
        for lam,mu,nu,a,b in allA(N):
            tot+=1
            if not cert('A',N,lam,mu,nu,a,b,routes): bad+=1
        fs=set(FAIL)
        print(N,routes,'allowed',tot,'uncert',bad,'distinct fail nodes',len(fs))
        if routes==('0','inf'): print('  sample fails',sorted(fs,key=lambda x:(x[1],str(x)))[:8])
