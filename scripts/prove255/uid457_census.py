# Reconcile Rick's "2225 shifted skew shapes 2<=n<=6" (brcheck.py) with Clio's 1243 (uid457).
# brcheck.py n enumerates (nu,lam) PAIRS of strict partitions with parts <= n, lam strictly inside nu.
import sys, itertools, collections
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake252/caseB')
from brpieri import boxes, se_rim
cum = 0
for n in range(2, 7):
    parts = [tuple(sorted(c, reverse=True)) for L in range(n + 1) for c in itertools.combinations(range(1, n + 1), L)]
    pairs = [(nu, lam) for nu in parts for lam in parts if boxes(lam) < boxes(nu)]
    S = {frozenset(boxes(nu) - boxes(lam)) for nu, lam in pairs}
    cum += len(pairs)
    rims = sum(se_rim(set(t)) == set(t) for t in S)
    print(f'n={n}: pairs {len(pairs)}, distinct box-sets {len(S)} (rims {rims}, non-rims {len(S)-rims}), cumulative pairs {cum}')
    if n == 6:
        print('  distinct box-sets with |theta| <= 13:', sum(len(t) <= 13 for t in S),
              '; non-rims by size:', sorted(collections.Counter(len(t) for t in S if se_rim(set(t)) != set(t)).items()))
