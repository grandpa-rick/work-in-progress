# Test: F_kappa(M) = C_kappa * w_M^2 * sum_{nu>=kappa} mu(kappa,nu) Psi_{nu,0}(M) ?
import sys, random
from fractions import Fraction as Fr
sys.path.insert(0,'../prove249c')
from eng import F, Psi, leq, subsets
n=int(sys.argv[1]); random.seed(int(sys.argv[2]) if len(sys.argv)>2 else 3)
w={i:Fr(random.randint(2,40),random.randint(2,40)) for i in range(1,n+1)}
S=[tuple(s) for s in subsets(n)]
le=lambda a,b: leq(list(a),list(b))
mu={}
for a in S:
    for b in sorted([b for b in S if le(a,b)],key=lambda b:(sum(b),len(b))):
        if a==b: mu[a,b]=1
        else: mu[a,b]=-sum(mu[a,c] for c in S if le(a,c) and le(c,b) and c!=b and (a,c) in mu)
def wp(M):
    p=Fr(1)
    for m in M: p*=w[m]
    return p
for k in S:
    rats=set()
    for M in S:
        lhs=F(list(k),list(M),w)
        rhs=wp(M)**2*sum(mu[k,v]*Psi(list(v),[0]*len(v),list(M),w) for v in S if le(k,v))
        if rhs==0 and lhs==0: continue
        rats.add(lhs/rhs if rhs!=0 else 'inf')
    print(k, 'ratio const' if len(rats)==1 else 'NOT const', list(rats)[:3] if len(rats)<=3 else len(rats))
