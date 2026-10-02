"""Reproduce MacBeth's per-brace table (D=Z/8 cyclic, H=Z/2, |G|=16): per L, (brace,D) pairs
with based == crit.  L given as set of units (lambda_d at a generator of D)."""
from collections import defaultdict
from braces import AbGroup, regular_subgroups, brace_ops
from criterion import run
for mods in [(16,), (2, 8)]:
    A, R, out = run(mods)
    tab = defaultdict(lambda: defaultdict(int))
    for r in out:
        D = r['D']
        if max(A.order(x) for x in D) != 8: continue
        g = next(x for x in D if A.order(x) == 8)
        co = {A.mul(k, g): k for k in range(8)}
        lam = R[r['brace']]
        L = tuple(sorted({co[lam[d][g]] for d in D}))
        Dlabel = tuple(sorted(D))[:3]
        tab[(L, Dlabel)][(r['based'], r['crit'])] += 1
    for (L, Dl), d in sorted(tab.items()):
        tot = sum(d.values()); ag = d[(True, True)] + d[(False, False)]
        print(f"Z{mods} D~{Dl} L={L}: pairs {tot}, based&crit {d[(True,True)]}, neither {d[(False,False)]}, "
              f"crit-but-unbased {d[(False,True)]}, agree {ag}/{tot}")
