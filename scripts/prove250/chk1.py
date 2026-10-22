# PROVE 250: check (i) Lemma N: s(N(sL)) = N(L); (ii) left action on gamma_nu = Psi_{nu,0}; (iii) support of ell
import sys, itertools
sys.path.insert(0,'../prove249c')
from fractions import Fraction as Fr
from eng import Psi, F, cc, ell, leq, subsets
n=int(sys.argv[1])
P=[5,7,11,13,17,19,23][:n]
w0={m:Fr(P[m-1]) for m in range(1,n+1)}
def sw(s,w):  # s on parameters
    w=dict(w)
    if s<n: w[s],w[s+1]=w[s+1],w[s]
    else: w[1]=1/w[1]
    return w
def sL(s,L):
    L=set(L)
    if s<n:
        a,b=s in L,(s+1) in L
        L-= {s,s+1}
        if a: L.add(s+1)
        if b: L.add(s)
    else: L^={1}
    return sorted(L)
def ema(s,w): return w[s]/w[s+1] if s<n else w[1]**-2   # e^{-alpha}
def comp(L): return [i for i in range(1,n+1) if i not in L]
def gam(nu,M,w): return Psi(nu,[0]*len(nu),M,w)
def N(L,w):
    wl=Fr(1)
    for m in L: wl*=w[m]**2
    return cc(L,w)/(wl*F(L,L,w)*gam(comp(L),comp(L),w))
# (i)
bad=0
for s in range(1,n+1):
    for L in subsets(n):
        if N(sL(s,L),sw(s,w0))!=N(L,w0): bad+=1
print('LemmaN mismatches',bad)
# (ii) for nu' and s with nu = smaller: s.gam_{nu'} = e^a gam_{nu'} - (e^a-1) gam_nu
bad=0;cnt=0;badS=0
for nup in subsets(n):
    for s in range(1,n+1):
        if s<n:
            if not ((s+1) in nup and s not in nup): continue
            nu=sorted([s if x==s+1 else x for x in nup])
        else:
            if 1 not in nup: continue
            nu=[x for x in nup if x!=1]
        ea=1/ema(s,w0)
        for M in subsets(n):
            lhs=gam(nup,sL(s,M),sw(s,w0))
            rhs=ea*gam(nup,M,w0)-(ea-1)*gam(nu,M,w0)
            cnt+=1
            if lhs!=rhs: bad+=1
print('left action gamma: instances',cnt,'mismatches',bad)
# (iii) support of ell: ell(kap,L)=0 unless L<=kap
bad=0
for k in subsets(n):
    for L in subsets(n):
        if not leq(L,k) and ell(k,L,n,w0)!=0: bad+=1
print('ell off-support nonzero',bad)
