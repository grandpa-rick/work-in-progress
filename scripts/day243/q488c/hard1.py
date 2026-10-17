import itertools,time
from resq0 import residue
t=time.time()
lam=(5,);mu=(5,);nu=(1,3,4,6);N=6
for a,b in itertools.combinations(range(1,7),2):
    r=residue(lam,mu,nu,a,b,N,seed=3)
    print((a,b),'Res!=0' if r!=0 else 'Res=0', str(r)[:150],flush=True)
print(time.time()-t)
