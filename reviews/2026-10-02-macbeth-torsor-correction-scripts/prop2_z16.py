"""MacBeth Prop. 2 is stated for D=Z/16 with lambda_g=5, 'L={1,5}'. Check: on Z/32 and Z/2xZ/16,
D cyclic of order 16, tabulate L (units at a generator), b parity, and criterion."""
import sys
from collections import defaultdict
from braces import AbGroup, regular_subgroups, brace_ops
from criterion import run
for mods in [tuple(map(int, a.split('x'))) for a in sys.argv[1:]]:
    A, R, out = run(mods)
    tab = defaultdict(lambda: defaultdict(int))
    for r in out:
        D = r['D']
        if max(A.order(x) for x in D) != 16: continue
        g = next(x for x in D if A.order(x) == 16)
        co = {A.mul(k, g): k for k in range(16)}
        lam = R[r['brace']]
        L = tuple(sorted({co[lam[d][g]] for d in D}))
        tab[L][(r['based'], r['crit'])] += 1
    print(f"Z{mods}: {len(R)} braces")
    for L, d in sorted(tab.items()):
        print(f"  D=Z/16 L={L}: (based,crit) {dict(d)}")
