import sys,importlib,ast; sys.path.insert(0,"../prove247b"); sys.path.insert(0,'../prove247')
import mp
from mp import Sys,subsets,EB,EC
mod=importlib.import_module(sys.argv[1])
N=4; fc=0; tot=0
for line in open('../prove247b/z4.txt'):
    c=int(line.split()[3]); zs=set(ast.literal_eval(line[line.index('['):line.rindex(']')+1]))
    for lam,mu,nu,a,b in mod.allA(N):
        if (a,b,lam,mu,nu) in zs or (a,b,mu,lam,nu) in zs:
            tot+=1; fc+=mod.start('A',N,lam,mu,nu,a,b,('eq',c))
print(sys.argv[1],'A N=4 zeros',tot,'falsecert',fc,flush=True)
for n in (3,4):
  s=[Sys(n,sd) for sd in (1,2)]
  for c in range(-1,2*n+1):
    for lam in subsets(n):
      for mu in subsets(n):
        for k in subsets(n):
          if EC(n,lam,mu,k) and all(mp.C(t,lam,mu,k,c)==0 for t in s):
            tot+=1; fc+=mod.start('C',n,lam,mu,k,None,None,('eq',c))
          for a in range(1,n+1):
            if EB(n,a,lam,mu,k) and all(mp.B(t,a,lam,mu,k,c)==0 for t in s):
              tot+=1; fc+=mod.start('B',n,lam,mu,k,a,None,('eq',c))
print(sys.argv[1],'total zeros',tot,'falsecert',fc,'PSTAT',mod.PSTAT)
