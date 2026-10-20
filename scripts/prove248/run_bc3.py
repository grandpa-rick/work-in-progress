import sys; sys.path.insert(0,"../prove247b"); from scert3 import *
import sys,time
for n in range(1,int(sys.argv[1])+1):
    t=time.time(); bc=tc=bb=tb=0; ex=[]
    for lam in subsets(n):
        for mu in subsets(n):
            if mu<lam: continue
            for k in subsets(n):
                if EC(n,lam,mu,k):
                    tc+=1
                    if not start('C',n,lam,mu,k,None,None): bc+=1; ex.append(('C',lam,mu,k))
                for a in range(1,n+1):
                    if EB(n,a,lam,mu,k):
                        tb+=1
                        if not start('B',n,lam,mu,k,a,None): bb+=1; ex.append(('B',lam,mu,k,a))
    print(n,'C allowed',tc,'uncert',bc,'| B allowed',tb,'uncert',bb,ex[:5],round(time.time()-t),'s',flush=True)
