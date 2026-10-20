import sys,importlib,time; sys.path.insert(0,"../prove247b"); sys.path.insert(0,'../prove247')
import mp
from mp import Sys,subsets,EC,EB
n=int(sys.argv[1]); mods=[importlib.import_module(m) for m in ('scert7','scert6m','scert5m2')]
s=[Sys(n,sd) for sd in (1,2)]
T=sys.argv[2] if len(sys.argv)>2 else 'C'
for c in range(-1,2*n+1):
    t=time.time(); fc=[0]*3; nz=0
    for lam in subsets(n):
      for mu in subsets(n):
        if mu<lam: continue
        for k in subsets(n):
          if T=='C' and EC(n,lam,mu,k) and all(mp.C(x,lam,mu,k,c)==0 for x in s):
            nz+=1
            for i,m in enumerate(mods): fc[i]+=m.start('C',n,lam,mu,k,None,None,('eq',c))
    print(n,T,'c',c,'zeros',nz,'falsecert [v3, P2-sign-mutant, P1p-mutant]',fc,[m.PSTAT for m in mods],round(time.time()-t),'s',flush=True)
