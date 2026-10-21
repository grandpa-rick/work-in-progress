# test: P^k_{lm}(c) = sum_nu c^nu_{lm}(1) P^k_{nu}(c); signs of c^nu(1) and of P^k_nu(U+1+t)
import sys, re, itertools, sympy as sp
from eng import Cfun, leq
from fractions import Fraction as Fr
n=int(sys.argv[1]); U=2*n-2; c=sp.Symbol('c'); t=sp.Symbol('t')
Ps={}
for fn in sys.argv[2:]:
    for line in open(fn):
        m=re.match(r'(\[.*?\]) (\[.*?\]) (\[.*?\]) d (\d+) deg (-?\d+).* P (.*)$',line.strip())
        lam,mu,kap=[tuple(eval(m.group(i))) for i in (1,2,3)]
        Ps[(lam,mu,kap)]=(int(m.group(4)),sp.sympify(m.group(6)))
def lim_c0(lam,mu,nu):
    # w->1 limit of C(lam,mu,nu;0) (raw structure constant), via Richardson
    a=[Fr(k*k+1,7) for k in range(1,n+1)]; vs=[]
    for e in (Fr(1,10**20),Fr(2,10**20),Fr(4,10**20)):
        w={i:1+a[i-1]*e for i in range(1,n+1)}; vs.append(Cfun(list(lam),list(mu),list(nu),0,n,w))
    return ((8*vs[0]-6*vs[1]+vs[2])/3).limit_denominator(10**6)
Uall=list(range(1,n+1)); subs=[tuple(S) for r in range(n+1) for S in itertools.combinations(Uall,r)]
sgnRich={}
for nu in subs:
  for kap in subs:
    if (( ),nu,kap) in Ps and Ps[((),nu,kap)][0]==0:
        Q=[x for x in sp.Poly(sp.expand(Ps[((),nu,kap)][1].subs(c,U+1+t)),t).all_coeffs() if x!=0]
        sgnRich[(nu,kap)]=sp.sign(Q[0])
# sign rule guess for Richardson: (-1)^{?}
from collections import Counter
cnt=Counter((sum(kap)-sum(nu))%2==0 and s==1 or (sum(kap)-sum(nu))%2==1 and s==-1 for (nu,kap),s in sgnRich.items()); print('Richardson sign = (-1)^{|kap|-|nu|}?',cnt)
cnt2=Counter(s for s in sgnRich.values()); print('Richardson signs',cnt2)
# structure constants and decomposition test
bad=0; tot=0; signs=Counter()
for lam in subs:
  for mu in subs:
    if lam>mu or lam==(): continue
    cs={nu:lim_c0(lam,mu,nu) for nu in subs if any(leq(M,nu) and leq(lam,M) and leq(mu,M) for M in subs)}
    for nu,v in cs.items():
        if v!=0: signs[(sp.sign(v)*(-1)**((sum(nu)-sum(lam)-sum(mu))%2))]+=1
    for kap in subs:
        if (lam,mu,kap) not in Ps or Ps[(lam,mu,kap)][0]!=0: continue
        tot+=1
        S=sum(sp.Rational(v.numerator,v.denominator)*Ps[((),nu,kap)][1] for nu,v in cs.items() if v!=0 and ((),nu,kap) in Ps and Ps[((),nu,kap)][0]==0)
        r=sp.simplify(Ps[(lam,mu,kap)][1]/S) if S!=0 else None
        if r is None or r.has(c): bad+=1
print('struct-const signs (x (-1)^{|nu|-|lam|-|mu|}):',signs)
print('decomposition proportional (const ratio) in',tot-bad,'of',tot)
print('exceptions:')
for lam in subs:
  for mu in subs:
    if lam>mu or lam==(): continue
    for nu in subs:
      if not any(leq(M,nu) and leq(lam,M) and leq(mu,M) for M in subs): continue
      v=lim_c0(lam,mu,nu)
      if v!=0 and sp.sign(v)*(-1)**((sum(nu)-sum(lam)-sum(mu))%2)==-1: print(lam,mu,nu,v)
      if v!=0 and sum(nu)<sum(lam)+sum(mu): print('  sub-degree nonzero:',lam,mu,nu,v)
