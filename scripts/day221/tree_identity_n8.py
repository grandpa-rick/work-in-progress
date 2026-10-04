# Day 221: history (Thm W weights) == telescoped tree form (Lemma 2) == Thm G, exact rationals, n<=8, l<=6, 4 t-values
from tree_identity_lib import *
def q(m,u): return sum(u**i for i in range(m))
def W(k,J,t):
    n=k+sum(J); r=Fr((-1)**len(J))*q(n,t)/q(k,t)
    for j in J: r*=q(k,t**j)
    return r
def parts(n,mx=None):
    if mx is None: mx=n
    if n==0: yield (); return
    for a in range(min(n,mx),0,-1):
        for r in parts(n-a,a): yield (a,)+r
cnt=0
for t in (Fr(3,5),Fr(7,3),Fr(-2),Fr(2)):
  for n in range(2,9):
    for lam in parts(n):
        l=len(lam)
        if l<2 or l>6: continue
        H=0; Tel=0
        for par in inc_trees(l):
            S=subtree(par,l); w=1; e=1
            for i in range(l):
                ch=[c for c in range(1,l) if par[c]==i]
                w*=W(lam[i],[sum(lam[j] for j in S[c]) for c in ch],t)
            for c in range(1,l): e*=t**(lam[par[c]]*sum(lam[j] for j in S[c]))-1
            H+=w; Tel+=e
        K=0
        for G in conn_graphs(l):
            w=1
            for a,b in G: w*=t**(lam[a]*lam[b])-1
            K+=w
        pre=Fr(1-t**n)
        for a in lam: pre/=(1-t**a)
        assert H==pre*Tel==pre*K,(lam,t); cnt+=1
print("history == telescoped tree == Thm G (t in 3/5,7/3,-2,2; n<=8, 2<=l<=6):",cnt)
