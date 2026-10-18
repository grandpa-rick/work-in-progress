import sys, itertools, time; sys.path.insert(0,'.')
from fast import *
N=int(sys.argv[1]); pairs=[tuple(map(int,p.split(','))) for p in sys.argv[2:]] or list(itertools.combinations(range(1,N+1),2))
Ss=subsets(N)
tri=[(l,m,nu) for i,l in enumerate(Ss) for m in Ss[i:] for nu in Ss if geq(nu,l) and geq(nu,m) and len(nu)-max(len(l),len(m))>=2]
def allowed(l,m,nu,a,b):
    rest=[c for c in range(1,N+1) if c not in (a,b)]
    for k in range(len(rest)+1):
        for Mpp in itertools.combinations(rest,k):
            M=tuple(sorted(Mpp+(a,b)))
            if geq(nu,M) and geq(Mpp,l) and geq(Mpp,m): return True
    return False
t0=time.time()
for a,b in pairs:
    B={(wt,sd):Basis(N,a,b,sd,wt) for wt in (0,1) for sd in (3,11)}
    st={}
    bad=[]
    for l,m,nu in tri:
        al=allowed(l,m,nu,a,b)
        z0=all(B[0,sd].resid(l,m,nu)==0 for sd in (3,11))
        zi=all(B[1,sd].resid(l,m,nu)==0 for sd in (3,11))
        k=(al,z0,zi); st[k]=st.get(k,0)+1
        if al==zi or (not al and not z0): bad.append((l,m,nu,'allowed' if al else 'notallowed','z0' if z0 else 'nz0','zinf' if zi else 'nzinf'))
    print('pair',a,b,'(allowed,Z0,Zinf):count',st,'%.0fs'%(time.time()-t0),flush=True)
    for x in bad[:20]: print('   ',x)
