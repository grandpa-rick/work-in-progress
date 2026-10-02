"""Pool (brace, D) pairs over all abelian (G,+) of a given order by triplet class;
compare class-level basepoint (exists member with b in 2D) with MacBeth's criterion
(exists member and section with delta_.(d) in Hom(Z/2,D) for all d)."""
import sys
from collections import defaultdict
from criterion import run
groups = {8: [(8,), (2, 4), (2, 2, 2)], 16: [(16,), (2, 8), (4, 4), (2, 2, 4)], 32: [(32,)]}
for order in map(int, sys.argv[1:]):
    cls = defaultdict(list)
    for mods in groups[order]:
        A, R, out = run(mods)
        print(f"  (G,+)=Z{mods}: {len(R)} braces, {len(out)} (brace,D) pairs", flush=True)
        for r in out:
            cls[r['inv']].append((mods, r))
    tab = defaultdict(lambda: defaultdict(int))
    bad = []
    for key, mem in cls.items():
        based = any(r['based'] for _, r in mem); crit = any(r['crit'] for _, r in mem)
        Lsize = mem[0][1]['Lsize']
        tab[(key[0], Lsize)][(based, crit)] += 1
        if based != crit: bad.append((key[0], Lsize, sorted({m for m, _ in mem}), len(mem), based, crit))
    print(f"=== |G|={order}: {len(cls)} triplet classes")
    for k, d in sorted(tab.items()):
        tot = sum(d.values()); agree = d[(True, True)] + d[(False, False)]
        print(f"  D={k[0]} |L|={k[1]}: classes {tot}, (based,crit) {dict(d)}, agree {agree}/{tot}")
    for b in bad: print("  MISMATCH class:", b)
