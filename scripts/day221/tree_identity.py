# Day 221 PROVE: check (i) increasing-tree identity with GENERIC edge weights x_ij;
# (ii) history sum with Thm W weights == telescoped tree formula == Conj G, symbolic t.
import itertools, sympy as sp, random
from fractions import Fraction as Fr
def inc_trees(l):  # parent arrays: par[c] < c for c=1..l-1 (0-indexed, root 0)
    for par in itertools.product(*[range(c) for c in range(1,l)]):
        yield (None,)+par
def subtree(par,l):
    S=[{i} for i in range(l)]
    for c in range(l-1,0,-1): S[par[c]] |= S[c]
    return S
def conn_graphs(l):
    E=list(itertools.combinations(range(l),2))
    for m in range(1<<len(E)):
        H=[E[i] for i in range(len(E)) if m>>i&1]
        p=list(range(l))
        def f(a):
            while p[a]!=a: a=p[a]
            return a
        for a,b in H: p[f(a)]=f(b)
        if len({f(a) for a in range(l)})==1: yield H
ok=0
for l in range(2,7):
    for trial in range(3):
        x={(i,j):Fr(random.randint(-9,9),random.randint(1,9)) for i in range(l) for j in range(i+1,l)}
        X=lambda i,j: x[(min(i,j),max(i,j))]
        L=0
        for par in inc_trees(l):
            S=subtree(par,l); w=1
            for c in range(1,l):
                pr=1
                for j in S[c]: pr*=X(par[c],j)
                w*=pr-1
            L+=w
        R=0
        for H in conn_graphs(l):
            w=1
            for e in H: w*=X(*e)-1
            R+=w
        assert L==R,(l,L,R); ok+=1
print("generic-x tree identity OK",ok)
t=sp.symbols('t')
q=lambda m,u=t: sum(u**i for i in range(m))
def W(k,J):
    n=k+sum(J); r=(-1)**len(J)*q(n)/q(k)
    for j in J: r*=q(k,t**j)
    return r
def parts(n,mx=None):
    if mx is None: mx=n
    if n==0: yield (); return
    for a in range(min(n,mx),0,-1):
        for r in parts(n-a,a): yield (a,)+r
cnt=0
for n in range(2,8):
    for lam in parts(n):
        l=len(lam)
        if l<2: continue
        H=0
        for par in inc_trees(l):
            S=subtree(par,l); w=1
            for i in range(l):
                ch=[c for c in range(1,l) if par[c]==i]
                w*=W(lam[i],[sum(lam[j] for j in S[c]) for c in ch])
            H+=w
        K=0
        for G in conn_graphs(l):
            w=1
            for a,b in G: w*=t**(lam[a]*lam[b])-1
            K+=w
        Gf=(1-t**n)*K
        for a in lam: Gf/=(1-t**a)
        assert sp.simplify(H-Gf)==0,lam; cnt+=1
print("history sum == Conj G, all lam |-n<=7:",cnt)
