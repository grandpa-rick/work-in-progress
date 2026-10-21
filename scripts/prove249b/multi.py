import sys,collections; sys.path.insert(0,"../prove248"); sys.path.insert(0,"../prove247b"); sys.path.insert(0,'../prove247')
import scert7 as S
MULTI=[]
orig=S.certify
def cert(Sx,n,a,b,cr,depth=0):
    Sc=S.clean(Sx,n,a,b)
    if len(Sc)>1: MULTI.append((n,a,b,cr,dict(Sc)))
    return orig(Sx,n,a,b,cr,depth)
S.certify=cert
N=int(sys.argv[1])
for x in S.allA(N): S.start('A',N,*x)
print(len(MULTI))
seen=set()
for n,a,b,cr,Sc in MULTI:
    keys=sorted(Sc)
    # describe what differs
    sig=tuple((k[1],k[2],k[3],k[4],k[5],k[6],v) for k,v in sorted(Sc.items()))
    if sig in seen: continue
    seen.add(sig)
    print('n',n,'a,b',a,b,'cr',cr)
    for k,v in sorted(Sc.items()):
        T,lam,th,mu,thp,nu,psi,d,e=k
        print('   ',v,T,'lam',lam,th,'mu',mu,thp,'dual',S.comp(n,nu),psi,'d',d,'e',e[1:])
