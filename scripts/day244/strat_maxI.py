# Wake 244: stratify N=5 (c)-failures by Clio's stratum max I <= b (I in {lam,mu}).
# Interpretation (Rick's, flagged to Clio): for triple (lam,mu,nu) and predicted pair (a,b)
# [pairs from n5scan.py: (a,b) in some k<=nu with k\{a,b} >= lam,mu], "in stratum" iff max(lam u mu) <= b.
import itertools, sys, ast
sys.argv=[sys.argv[0],'1']
from engine import subsets, geq
from coverage import cover
N=5
S=subsets(N)
def rm(M,a,b): return tuple(x for x in M if x not in (a,b))
def pairs(l,m,nu):
    P=[]
    for k in S:
        if not geq(nu,k): continue
        for a,b in itertools.combinations(k,2):
            r=rm(k,a,b)
            if geq(r,l) and geq(r,m) and (a,b) not in P: P.append((a,b))
    return P
fails=ast.literal_eval(open('n5fails_partial.txt').read())
tot,ct,cs,unc=cover(N)
allt=[]
for i,l in enumerate(S):
    for m in S[i:]:
        for nu in S:
            if geq(nu,l) and geq(nu,m) and len(nu)-max(len(l),len(m))>=2: allt.append((l,m,nu))
def cls(t):
    l,m,nu=t; mx=max(l+m) if (l+m) else 0
    P=pairs(*t); ins=[mx<=b for a,b in P]
    return any(ins), all(ins)
def tab(name,T):
    c=[cls(t) for t in T]
    print(f"{name}: n={len(T)}  some-pair-in-stratum={sum(a for a,b in c)}  all-pairs-in-stratum={sum(b for a,b in c)}  no-pair-in-stratum={sum(not a for a,b in c)}")
tab('all jump>=2',allt); tab('uncovered by Thm9/9prime',unc)
tab('c-fails',fails); tab('uncovered non-fails',[t for t in unc if t not in fails])
for t in fails:
    l,m,nu=t; print(t,'max I=',max(l+m),'pairs',pairs(*t))
