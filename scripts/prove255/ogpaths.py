# PROVE 255 OG twin: type B_n, spin weight omega_n=(1/2,..,1/2); covers/bonds from root data; LS count vs og.chi_dem and Omega(P,j+1)
import sys, itertools
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0,'../wake253'); sys.path.insert(0,'../wake255')
import og
from models_from_reciprocity import count
import lspaths as L
def rootsB(n):
    R=[]
    for i in range(n):
        for k in range(i+1,n):
            for s in (1,-1):
                v=[0]*n; v[i]=1; v[k]=s; R.append((tuple(v),tuple(v)))
        v=[0]*n; v[i]=1; c=[0]*n; c[i]=2; R.append((tuple(v),tuple(c)))
    return R
def setupB(n):
    R=rootsB(n); h=Fr(1,2); V=list(itertools.product([h,-h],repeat=n))
    D={v:L.dim(v,R) for v in V}; covers=defaultdict(list)
    for v in V:
        for (a,c) in R:
            u=L.refl(v,a,c)
            if D[u]==D[v]+1: covers[v].append((u,abs(L.pair(v,c))))
    up={v:{v} for v in V}
    for v in sorted(V,key=lambda x:-D[x]):
        for (u,b) in covers[v]: up[v]|=up[u]
    return R,V,D,covers,up
N=int(sys.argv[1]); J=int(sys.argv[2]); S=setupB(N)
print('bonds',sorted({b for c in S[3].values() for (u,b) in c}))
# check dim formula matches parts convention
tot=fails=0
for Lr in range(N+1):
    for c in itertools.combinations(range(1,N+1),Lr):
        kap=tuple(sorted(c,reverse=True)); h=Fr(1,2)
        kv=tuple(-h if (N-i) in c else h for i in range(N))
        assert S[2][kv]==sum(c)
        for j in range(1,J+1):
            P=L.lspaths(kv,N,j,S); Om=count(kap,range(0,j+1),None,False); H=og.chi_dem(kap,N,j)
            tot+=1; ok=(len(P)==Om==H); fails+=not ok
            if not ok: print('FAIL',kap,j,len(P),Om,H)
print(f'OG n={N} j<={J}: {tot} cases, failures {fails}')
