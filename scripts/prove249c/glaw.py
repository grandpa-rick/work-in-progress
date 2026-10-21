import re, sys
from collections import defaultdict
tab=defaultdict(int)
for fn in sys.argv[1:]:
  for line in open(fn):
    m=re.match(r'\((\d+), (\d+)\) (\(.*?\)) (\(.*?\)) (\(.*?\)) d (\d+) deg (-?\d+)',line)
    if not m: continue
    tl=lambda x: x if isinstance(x,tuple) else (x,)
    lam,mu,nu=[tl(eval(m.group(i))) for i in (3,4,5)]
    tab[(int(m.group(6))-int(m.group(7)))-(sum(lam)+sum(mu)-sum(nu))]+=1
print(dict(tab))
