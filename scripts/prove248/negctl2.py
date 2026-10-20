import sys; sys.path.insert(0,"../prove247b"); from scert2 import *
import ast,re,sys
for fn,N in (('../prove247b/z4.txt',4),('../prove247b/z5.txt',5)):
    for line in open(fn):
        c=int(line.split()[3]); zs=ast.literal_eval(line[line.index('['):line.rindex(']')+1])
        zset=set(zs); falsecert=[]; ncert=0; nz=0
        for lam,mu,nu,a,b in allA(N):
            key=(a,b,lam,mu,nu); key2=(a,b,mu,lam,nu)
            isz = key in zset or key2 in zset
            ok=start('A',N,lam,mu,nu,a,b,('eq',c))
            if isz and ok: falsecert.append(key)
            if not isz: nz+=1; ncert+=ok
        print(N,'c',c,'zeros',len(zs),'FALSE-CERT',len(falsecert),falsecert[:3],'| nonzero certified',ncert,'/',nz,flush=True)
