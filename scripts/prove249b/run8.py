import sys,time,collections
import scert8 as S
S.MUT=sys.argv[3] if len(sys.argv)>3 else None
N=int(sys.argv[1]); mode=sys.argv[2]
MULTI=collections.Counter(); orig=S.certify
def cert(Sx,n,a,b,cr,depth=0):
    Sc=S.clean(Sx,n,a,b)
    if len(Sc)>1: MULTI['multi']+=1
    if Sc: MULTI['nodes']+=1
    return orig(Sx,n,a,b,cr,depth)
S.certify=cert
t=time.time()
if mode=='run':
    bad=[];tot=0
    for x in S.allA(N):
        tot+=1
        if not S.start('A',N,*x): bad.append(x)
    print('v8 A N',N,'allowed',tot,'uncert',len(bad),bad[:5],round(time.time()-t),'s',dict(MULTI),S.PSTAT,flush=True)
else:
    import ast
    fc=0;tot=0;nzc=0;nz=0
    for fn,NN in (('../prove247b/z4.txt',4),('../prove247b/z5.txt',5)):
        if NN>N: continue
        for line in open(fn):
            c=int(line.split()[3]); zs=set(ast.literal_eval(line[line.index('['):line.rindex(']')+1]))
            for lam,mu,nu,a,b in S.allA(NN):
                isz=(a,b,lam,mu,nu) in zs or (a,b,mu,lam,nu) in zs
                ok=S.start('A',NN,lam,mu,nu,a,b,('eq',c))
                if isz: tot+=1; fc+=ok
                else: nz+=1; nzc+=ok
    print('negctl v8 MUT',S.MUT,'zeros',tot,'FALSE-CERT',fc,'| nonzero certified',nzc,'/',nz,round(time.time()-t),'s',flush=True)
