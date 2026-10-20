import sys,importlib,time; sys.path.insert(0,"../prove247b"); sys.path.insert(0,'../prove247')
import mp
from mp import Sys,subsets,EC
n=6; mods=[importlib.import_module(m) for m in ('scert7','scert7m')]
s=[Sys(n,sd) for sd in (1,2)]
nus=[k for k in subsets(n) if not ({2,3,4,5}&set(k))]
print('nus',nus,flush=True)
for c in range(-1,2*n+1):
    t=time.time(); fc=[0,0]; nz=0; tot=0
    for lam in subsets(n):
      for mu in subsets(n):
        if mu<lam: continue
        for k in nus:
          if EC(n,lam,mu,k):
            tot+=1
            if all(mp.C(x,lam,mu,k,c)==0 for x in s):
              nz+=1
              for i,m in enumerate(mods): fc[i]+=m.start('C',n,lam,mu,k,None,None,('eq',c))
    print(n,'c',c,'inst',tot,'zeros',nz,'falsecert [v4, P2-mutant]',fc,[m.PSTAT.get('P2',0) for m in mods],round(time.time()-t),'s',flush=True)
