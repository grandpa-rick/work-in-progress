from cert import *
import collections
N=int(sys.argv[1])
FAIL.clear()
for lam,mu,nu,a,b in allA(N): cert('A',N,lam,mu,nu,a,b,('0','inf'))
fs=set(FAIL); c=collections.Counter((f[0],f[1]) for f in fs); print(c)
mins=[f for f in fs if f[1]==min(x[1] for x in fs)]
for f in sorted(mins,key=str)[:40]:
    T,n,lam,mu,nu,a,b=f; E=supp(T,n,lam,mu,nu,a,b)
    print(f, 'E',E)
