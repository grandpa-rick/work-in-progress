import sys,collections; sys.path.insert(0,"../prove247b"); sys.path.insert(0,'../prove247')
import scert6 as S
CNT=collections.Counter(); MULTI=collections.Counter()
orig=S.certify
def cert(Sx,n,a,b,cr,depth=0):
    Sc=S.clean(Sx,n,a,b)
    if n==0 and Sc: CNT[sum(Sc.values())]+=1
    if Sc: MULTI[min(len(Sc),5)]+=1
    return orig(Sx,n,a,b,cr,depth)
S.certify=cert
N=int(sys.argv[1]); T=sys.argv[2]
if T=='A':
    for x in S.allA(N): S.start('A',N,*x)
else:
    from mp import subsets,EC
    for lam in subsets(N):
        for mu in subsets(N):
            if mu<lam: continue
            for k in subsets(N):
                if EC(N,lam,mu,k): S.start('C',N,lam,mu,k,None,None)
print(T,N,'base signed counts',sorted(CNT.items())[:20]); print('node sizes',sorted(MULTI.items()))
