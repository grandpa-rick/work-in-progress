# detect psi_opts ties (>=2 options at extremal order) and classify
import sys,collections; sys.path.insert(0,"../prove248"); sys.path.insert(0,"../prove247b"); sys.path.insert(0,'../prove247')
import scert7 as S
TIES=collections.Counter()
orig=S.psi_opts
def blocks(lam):
    b=[];j=0;r=len(lam)
    while j<r:
        q=j
        while q+1<r and lam[q+1]==lam[q]+1: q+=1
        b.append((j,q)); j=q+1
    return b
def po(lam,th,n,lim):
    out=orig(lam,th,n,lim)
    if out:
        ords=[o[0] for o in out]; ext=min(ords) if lim==0 else max(ords)
        k=sum(1 for o in ords if o==ext)
        if k>1:
            # P1'/P2 outputs are intentionally multi; skip those (all same order by construction)
            nb=len(blocks(lam))
            tag='P1p/P2' if (S.LEMMAP and ((len(lam)==1 and lam[0]>=2 and th[0]==lam[0] and lim==0) or (nb==1 and len(set(th))==1 and lam[0]>=2 and th[0]==lam[0] and lim==0))) else 'TIE'
            TIES[(tag,lim,lam,th,nb)]+=1
    return out
S.psi_opts=po
T=sys.argv[1]; N=int(sys.argv[2])
if T=='A':
    for x in S.allA(N): S.start('A',N,*x)
else:
    from mp import subsets,EC,EB
    for lam in subsets(N):
        for mu in subsets(N):
            if mu<lam: continue
            for k in subsets(N):
                if T=='C' and EC(N,lam,mu,k): S.start('C',N,lam,mu,k,None,None)
agg=collections.Counter()
for (tag,lim,lam,th,nb),v in TIES.items(): agg[(tag,lim,nb)]+=v
print(T,N,sorted(agg.items()))
for k,v in sorted(TIES.items(),key=lambda z:-z[1]):
    if k[0]=='TIE': print(k,v)
