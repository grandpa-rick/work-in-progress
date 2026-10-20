import sys; sys.path.insert(0,'../prove247')
import mp
from mp import Sys,subsets,EB,EC
import sys; sys.path.insert(0,"../prove247b"); import scert5 as scert
from scert5 import start
for n in (2,3,4):
  s=[Sys(n,sd) for sd in (1,2)]
  for c in range(-1,2*n+1):
    fcC=fcB=0; zC=zB=0; cC=cB=0; ex=[]
    for lam in subsets(n):
      for mu in subsets(n):
        for k in subsets(n):
          if EC(n,lam,mu,k):
            z=all(mp.C(t,lam,mu,k,c)==0 for t in s); ok=start('C',n,lam,mu,k,None,None,('eq',c))
            zC+=z; fcC+= (z and ok); cC+=ok
            if z and ok: ex.append(('C',lam,mu,k))
          for a in range(1,n+1):
            if EB(n,a,lam,mu,k):
              z=all(mp.B(t,a,lam,mu,k,c)==0 for t in s); ok=start('B',n,lam,mu,k,a,None,('eq',c))
              zB+=z; fcB+=(z and ok); cB+=ok
              if z and ok: ex.append(('B',lam,mu,k,a))
    print(n,c,'C zeros',zC,'falsecert',fcC,'| B zeros',zB,'falsecert',fcB,ex[:4],flush=True)
