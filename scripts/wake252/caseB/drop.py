# (Ind) mechanism: supp[O_D] = {rho: theta=kappa/rho rim with <=2 boxes per row/col} (BR Cor 5.6, p=1, dualised; checked support.py)
# => any chain meets theta in <=2 boxes => r'(rho) >= r'(kappa)-2.  Test: which theta realise drop -2?
import sys, collections
from support import parse, bx, rp
hits = collections.Counter()
for m in sys.argv[1:]:
    for k, K in parse(f'out_gate6_m{m}.txt'):
        if len(k) < 2: continue
        for rho in K:
            th2 = {(a, b) for (a, b) in bx(k) - bx(rho) if a >= 1}
            hdom = any((a, b + 1) in th2 for (a, b) in th2); vdom = any((a + 1, b) in th2 for (a, b) in th2)
            hits[(rp(rho) - rp(k), 'hdom' if hdom else '-', 'vdom' if vdom else '-')] += 1
for key, v in sorted(hits.items()): print(key, v)
