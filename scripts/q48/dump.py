import sys
src=open('../q47/ls12_cor34_check.py').read().split('allok = True')[0]
sys.argv=['x',sys.argv[1]]
exec(src)
def part(c): return tuple(sorted([m-i for i in range(m) if c[i]==1],reverse=True))  # c = -v_lambda
from collections import defaultdict
bad=0; sizes=Counter()
for k in strict:
  z=minrep[tuple(-x for x in v_of(k))]
  for r in strict:
    o=ours(k,r)
    if not o or not o[0]: continue
    w=minrep[tuple(-x for x in v_of(r))]
    ps=[p for p in paths if cleq(p[0][0],coset[w]) and downp(w,p)==z]
    th=sorted(D(k)-D(r))
    print('kappa',k,'rho',r,'theta',th)
    for p in ps:
      print('   path',[part(c) for c in p[0]],'mono',tuple(-x for x in endpoint(p)))
    for mono,c in o[0].items(): sizes[c]+=1
print(sizes)
