import time
from sy import *
for (a,b,l,m,nu,c) in [(1,2,(),(3,),(1,2,4),0),(1,2,(),(2,),(1,3,4),1),(1,3,(),(2,),(1,3,4),1),(1,4,(),(2,),(1,3,4),1),
                       (1,2,(),(3,),(1,2,4),1),(1,2,(),(3,),(1,2,4),-1)]:
    t=time.time(); s=Sys(4,0,a,b); r=A(s,a,b,l,m,nu,c)
    print((a,b,l,m,nu,c),'residue == 0:', r==FLD.zero, '%.0fs'%(time.time()-t), flush=True)
