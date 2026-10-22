import sys, itertools, sympy as sp
sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/prove249d')
from dem import Di
from collections import defaultdict
n=int(sys.argv[1]); A=int(sys.argv[2]) if len(sys.argv)>2 else 4
for a in range(0,A+1):
    f={}
    for S in itertools.product((0,1),repeat=n):
        k=tuple(-a-s for s in S); f[k]=f.get(k,0)+1
    for i in range(n,0,-1): f=Di(f,i,n)
    G=defaultdict(int)
    for k,c in f.items(): G[k[1:]]+=c
    G={k:v for k,v in G.items() if v}
    # check symmetry and decompose: for symmetric G, coefficient of e^{-b w'} e_j  = coeff of monomial with exponent pattern (-b-1)^j (-b)^{n-1-j}
    keys=set(G)
    sym=all(G.get(tuple(sorted(k,reverse=True)),0)==v for k,v in G.items())
    # decomposition attempt: greedy by dominant monomials
    rem=dict(G); dec={}
    while rem:
        k=max(rem, key=lambda t:(sorted(t,reverse=True)))  # pick a term
        kk=tuple(sorted(k,reverse=True)); c=rem.get(kk,0)
        if c==0: kk=k; c=rem[k]
        # need kk of form (-b)^{n-1-j}(-b-1)^j
        vals=set(kk)
        if not (len(vals)<=2 and (len(vals)==1 or max(vals)-min(vals)==1)):
            dec['NONBASIS',kk]=c; del rem[k]; continue
        b=-max(vals); j=sum(1 for v in kk if v==-b-1)
        dec[(b,j)]=c
        for S in itertools.product((0,1),repeat=n-1):
            if sum(S)==j:
                t=tuple(-b-s for s in S); rem[t]=rem.get(t,0)-c
                if rem[t]==0: del rem[t]
    print('n',n,'a',a,'sym',sym, dict(sorted(dec.items(),key=str)))
