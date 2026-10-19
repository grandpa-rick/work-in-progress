import scert,sys
from scert import *
N=int(sys.argv[1])
for L in ((0,),(1,)):
    scert.LIMS=L; bad=0;tot=0
    for lam,mu,nu,a,b in allA(N):
        tot+=1; bad+= not start('A',N,lam,mu,nu,a,b)
    print(N,'policy',L,'uncert',bad,'/',tot,flush=True)
