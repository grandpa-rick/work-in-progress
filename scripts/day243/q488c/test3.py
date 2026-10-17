import sympy as sp,sys,time
from engine import subsets, geq
import engine2, engine3
N=int(sys.argv[1]);t=time.time()
for M in subsets(N):
  for I in subsets(N):
    if not geq(M,I): continue
    a=engine2.A01(N,M,I); b=engine3.A01(N,M,I)
    if sp.cancel(a[0]-b[0])!=0 or sp.cancel(a[1]-b[1])!=0: print('MISMATCH',M,I)
print('ok',N,time.time()-t)
