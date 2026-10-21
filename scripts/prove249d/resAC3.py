import sys, random
from fractions import Fraction as Fr
sys.path.insert(0,'../prove249c')
from eng import ell, subsets
for N in (4,5):
    random.seed(4); base={i:Fr(random.randint(2,40),random.randint(2,40)) for i in range(1,N+1)}
    a,b=N-1,N
    def R(nu,L):
        out=[]
        for dl in (Fr(1,10**25),Fr(2,10**25)):
            w=dict(base); w[b]=1/(w[a]*(1+dl)); h=1-1/(w[a]*w[b]); out.append(ell(nu,L,N,w)*h)
        return 2*out[0]-out[1]
    v=base[a]; G=(1-v**2)*(1-v**-2)
    for k in range(1,N-1): G*=(1-v*base[k])*(1-base[k]/v)
    wsm={k:base[k] for k in range(1,N-1)}
    nu=list(range(1,N+1)); L=[a,b]
    K=R(nu,L)/ell(list(range(1,N-1)),[],N-2,wsm)
    print(N,'K*v*G =',float(K*v*G))
