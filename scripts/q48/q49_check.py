import sys; sys.path.insert(0,'prove262')
from q46 import *
import sympy as sp
m=int(sys.argv[1]); 
for a in (list(range(m,0,-1)), [2*(m-i) for i in range(m)][:-1]+[1] if False else [2*i-1 for i in range(m,0,-1)]):
  bad=0
  for k in strict(m):
    R=[r for r in strict(m) if sub(r,k)]
    binv={r: sum(-(1 if (m-i) in r else -1)*a[i] for i in range(m)) for r in R}  # (eps_r,a), eps_r=-1 at pos m+1-k
    S=sorted({sum(e*x for e,x in zip(eps,a)) for eps in itertools.product((1,-1),repeat=m)})
    def rk(rs):
      if not rs: return 0
      M=sp.Matrix([[sp.Rational(*fspec(chw(r,j,m),a,7)) for r in rs] for j in range(len(S))]); return M.rank()
    bins=sorted(set(binv.values()))
    gr=[rk([r for r in R if binv[r]>=b])-rk([r for r in R if binv[r]>b]) for b in bins]
    if any(g!=1 for g in gr): bad+=1; print('FAIL',k,gr)
    # sign check: min exponent of f_r equals bin(r)
  print('m',m,'a',a,'kappas with some gr_b != 1:',bad)
