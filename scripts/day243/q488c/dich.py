"""For each jump>=2 triple at N: full-q residue on H_ab (numeric w, symbolic q): q-valuation at 0 and degree at infinity,
compared with e-1, e = dim + codim nu - codim lam - codim mu (dim = N(N+1)/2, codim = sum of elements)."""
import itertools, random, sys, sympy as sp, time
from engine import geq, subsets
from fullc import A, q, eps
N=int(sys.argv[1]); a,b=int(sys.argv[2]),int(sys.argv[3]); t0=time.time()
random.seed(5)
wv=[sp.Rational(random.randint(2,97),random.randint(2,97))+random.randint(0,3) for c in range(N)]
wv[b-1]=1/(wv[a-1]*(1-eps))
S=subsets(N); order=sorted(S,key=lambda T:(len(T),sum(T)))
Am={(M,I):A(N,M,I,wv) for M in S for I in S}
print('A built',time.time()-t0,flush=True)
dim=N*(N+1)//2
for i,l in enumerate(S):
  for m in S[i:]:
    c={}
    for M in order:
      v=Am[M,l]*Am[M,m]
      c[M]=sp.cancel((v-sum(c[k]*Am[M,k] for k in c if Am[M,k]!=0 and c[k]!=0))/Am[M,M])
    for nu in S:
      if not(geq(nu,l) and geq(nu,m)) or len(nu)-max(len(l),len(m))<2: continue
      r=sp.cancel(sp.cancel(c[nu]*eps).subs(eps,0))
      if r==0: print(l,m,nu,'NO POLE'); continue
      n_,d_=sp.fraction(r); Pn=sp.Poly(n_,q); Pd=sp.Poly(d_,q)
      v0=min(x[0] for x in Pn.monoms())-min(x[0] for x in Pd.monoms())
      dinf=Pn.degree()-Pd.degree()
      e=dim+sum(nu)-sum(l)-sum(m)
      print(l,m,nu,'val0',v0,'deg_inf',dinf,'e-1',e-1, 'Q0' if v0==1 else '', 'QINF' if dinf==e-1 else '',flush=True)
print('time',time.time()-t0)
