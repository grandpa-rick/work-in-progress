# check: (-1)^|k| 2^m H(-M/2) = sum_S Dbar_{w_max}(e^{lam'_S})(1), lam'_S=(a-i+s_i), each term >=0 when a>=m
import sys, itertools
sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/prove249d')
from dem import Di, word
from pq import polys, x
import sympy as sp
def w0A(n): return [j for i in range(1,n) for j in range(i,0,-1)]
def Dbar(f,i,n):
    g=dict(Di(f,i,n))
    for k,v in f.items(): g[k]=g.get(k,0)-v
    return {k:v for k,v in g.items() if v}
def atom(kap,n,lam):
    f={tuple(lam):1}
    for i in reversed(word(set(kap),n)+w0A(n)): f=Dbar(f,i,n)
    return sum(f.values())
n=int(sys.argv[1]); ok=tot=0
for r in range(1,n+1):
  for kap in itertools.combinations(range(1,n+1),r):
    S_,E,O,p,q=polys(kap)
    for a in range(n-2,n+3):
        terms={S:atom(kap,n,[a-(i+1)+S[i] for i in range(n)]) for S in itertools.product((0,1),repeat=n)}
        lhs=(-1)**sum(kap)*2**n*E.subs(x,-(2*a+1))
        good = lhs==sum(terms.values()); tot+=1; ok+=good
        neg=[v for v in terms.values() if v<0]
        if a>=n: assert all(v>=0 for v in terms.values()) and terms[tuple([1]*n)]>=1, (kap,a,terms)
        if not good or (a<n and neg): print(kap,'a',a,'match',good,'neg terms',neg[:5])
print('identity',ok,'/',tot)
