import scert, sys
from scert import *
orig=scert.certify; SEEN={}
def certify(S,n,a,b,cr,depth=0):
    r=orig(S,n,a,b,cr,depth)
    S2=scert.clean(S,n,a,b)
    if len(S2)>=2:
        k=(n,a,b,cr,tuple(sorted(S2.items())))
        SEEN[k]=r
    return r
scert.certify=certify
N=int(sys.argv[1])
for lam,mu,nu,a,b in allA(N): start('A',N,lam,mu,nu,a,b)
print(len(SEEN),'multi nodes; certified',sum(SEEN.values()))
for k,r in list(SEEN.items())[:12]:
    n,a,b,cr,items=k; print('n',n,'ab',a,b,cr,'cert',r)
    for key,v in items: print('    ',v,key[0],'lam',key[1],key[2],'mu',key[3],key[4],'nu',key[5],'psi',key[6],'d',key[7],'e',key[8][1:])
