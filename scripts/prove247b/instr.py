import scert, collections, sys
from scert import *
# instrument: record (n, lim chosen at success, #summands entering) along successful path
orig=scert.certify
LOG=collections.Counter(); MAXS=collections.Counter()
def certify(S,n,a,b,cr,depth=0):
    r=orig(S,n,a,b,cr,depth)
    if r: MAXS[len(S)]+=1
    return r
scert.certify=certify
N=int(sys.argv[1])
for lam,mu,nu,a,b in allA(N): start('A',N,lam,mu,nu,a,b)
print('summand-count histogram over successful certify calls',sorted(MAXS.items()))
