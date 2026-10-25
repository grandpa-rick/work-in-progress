# Independent brute-force check (all permutations) that the only unmarked-isomorphic pair of shifted-diagram posets
# with |kappa|<=7, parts<=6, is (3)~(2,1).  Canonical form = lexicographically least relation matrix over all orderings.
import itertools, sys
sys.argv = [sys.argv[0], '6']; sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake253')
from og import strict_parts, boxes, nm
def le(a, b): return a[0] <= b[0] and a[1] <= b[1]
def canon(k):
    B = sorted(boxes(k)); return min(tuple(le(p[i], p[j]) for i in range(len(B)) for j in range(len(B))) for p in itertools.permutations(B))
cl = {}
for k in strict_parts(6):
    if k and sum(k) <= 7: cl.setdefault((sum(k), canon(k)), []).append(nm(k))
print('classes with >1 member:', [v for v in cl.values() if len(v) > 1], ' total shapes', sum(map(len, cl.values())))
