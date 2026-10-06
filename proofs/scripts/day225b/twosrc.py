from nonsym import *
from test1pt import comps
import sys
a = int(sys.argv[1])
for k in range(1, a):
  for x in range(1, 5):
    for y in range(1, 5):
        n = x + y; d = n - a
        if d < 0: continue
        rows = []
        for beta in comps(d, a):
            gam = [1 + b for b in beta]; gam[0] -= x; gam[k] -= y
            v = ctmono(tuple(gam))
            if v != 0: rows.append((beta, sp.factor(v)))
        print(f'a={a} src2=pos{k+1} x={x} y={y}:', rows)
