# Are {H_rho : rho strictly inside kappa} linearly independent? If yes for kappa, (star)+(hyperplane) would force Q30 at kappa.
import sys, itertools, sympy as sp
sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/prove249d')
from dem import chi as chiC
x=sp.Symbol('x'); m=int(sys.argv[1])
parts=[tuple(sorted(c,reverse=True)) for L in range(m+1) for c in itertools.combinations(range(1,m+1),L)]
H={k:sp.Poly(sp.interpolate([(j,chiC(set(k),m,j)[0]) for j in range(sum(k)+1)],x) if sum(k) else sp.Integer(1),x) for k in parts}
sub=lambda r,k: all(i<len(k) and r[i]<=k[i] for i in range(len(r))) and r!=k
dep=[]
for k in parts:
    R=[r for r in parts if sub(r,k)]
    if not R: continue
    D=sum(k)
    M=sp.Matrix([[H[r].coeff_monomial(x**i) for i in range(D)] for r in R])
    if M.rank()<len(R): dep.append(k)
print('m',m,'kappa with dependent {H_rho: rho<kappa}:',dep)
