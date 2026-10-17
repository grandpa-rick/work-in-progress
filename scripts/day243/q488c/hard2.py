import itertools,time,sys
from resnum import residue
t=time.time()
lam=eval(sys.argv[1]);mu=eval(sys.argv[2]);nu=eval(sys.argv[3]);N=int(sys.argv[4])
for a,b in itertools.combinations(range(1,N+1),2):
    r=residue(lam,mu,nu,a,b,N,seed=3)
    print((a,b),'POLE' if r!=0 else '-', flush=True)
print('time',time.time()-t)
