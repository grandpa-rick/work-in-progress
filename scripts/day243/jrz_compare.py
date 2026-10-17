# Day 243 (2026-10-17): JRZ 0908.0083 quantum quasi-shuffle product, implemented from the printed
# inductive formula (JRZ Prop 3, eq. (1)) with V=A_{phi,+}, sigma(u_a(x)u_b)=q^{(a,b)}u_b(x)u_a, m(u_a(x)u_b)=q^{phi(a,b)}u_{a+b},
# compared to Rick Day 241 QSym_{q,phi} weight formula. Exponents of q tracked exactly. Run with arg ctl1/ctl2 for controls.
import itertools, random, sys
CTL=sys.argv[1] if len(sys.argv)>1 else ''
from functools import lru_cache
from collections import defaultdict
B=[[2,-1],[-1,2]]  # symmetric form (sl3)
random.seed(1)
PHI=[[random.randint(-3,3) for _ in range(2)] for _ in range(2)]  # arbitrary bilinear, nonsymmetric
def form(a,b): return sum(a[i]*B[i][j]*b[j] for i in range(2) for j in range(2))
def phi(a,b): return sum(a[i]*PHI[i][j]*b[j] for i in range(2) for j in range(2))
def add(d,w,e,c=1):
    d[(w,e)]+=c
@lru_cache(None)
def jrz(U,V):
    out=defaultdict(int)
    if not U: out[(V,0)]+=1; return frozenset(out.items())
    if not V: out[(U,0)]+=1; return frozenset(out.items())
    for (w,e),c in jrz(U,V[:-1]): out[(w+(V[-1],),e)]+=c
    u=U[-1]; s=sum(form(u,v) for v in V)+(1 if CTL=='ctl2' else 0)
    for (w,e),c in jrz(U[:-1],V): out[(w+(u,),e+s)]+=c
    s2=sum(form(u,v) for v in V[:-1])+(phi(V[-1],u) if CTL=='ctl1' else phi(u,V[-1])); m=tuple(x+y for x,y in zip(u,V[-1]))
    for (w,e),c in jrz(U[:-1],V[:-1]): out[(w+(m,),e+s2)]+=c
    return frozenset((k,c) for k,c in out.items() if c)
def rick(A,Bw):
    out=defaultdict(int); r,s=len(A),len(Bw)
    for mlen in range(max(r,s),r+s+1):
        for F in itertools.combinations(range(mlen),r):
            for G in itertools.combinations(range(mlen),s):
                if set(F)|set(G)!=set(range(mlen)): continue
                e=0; word=[[0,0] for _ in range(mlen)]
                for i,f in enumerate(F): word[f]=[x+y for x,y in zip(word[f],A[i])]
                for j,g in enumerate(G): word[g]=[x+y for x,y in zip(word[g],Bw[j])]
                for i,f in enumerate(F):
                    for j,g in enumerate(G):
                        if g<f: e+=form(Bw[j],A[i])
                        if g==f: e+=phi(A[i],Bw[j])
                out[(tuple(map(tuple,word)),e)]+=1
    return frozenset((k,c) for k,c in out.items() if c)
letters=[(1,0),(0,1),(1,1),(2,1)]
words=[()]+[w for L in (1,2,3) for w in itertools.product(letters,repeat=L)]
ok=bad=0
for U in words:
    for V in words:
        if len(U)+len(V)>5: continue
        if jrz(U,V)==rick(U,V): ok+=1
        else: bad+=1
print(CTL or "main","PHI",PHI,"agree",ok,"disagree",bad)
