# PROVE 252 Lemma (chain drop), pure combinatorics, no Pieri input: for all strict rho <= kappa (parts in [n]) with
# theta = D(kappa)\D(rho) having no strictly-SE pair and <=2 boxes per row and column:  r'(rho) >= r'(kappa) - 2.
import sys, itertools, collections
from support import rp, bx
for n in range(1, int(sys.argv[1]) + 1):
    S = [tuple(sorted(c, reverse=True)) for L in range(n + 1) for c in itertools.combinations(range(1, n + 1), L)]
    tot = bad = 0
    for k in S:
        for r in S:
            if r == k or len(r) > len(k) or any(r[i] > k[i] for i in range(len(r))): continue
            th = bx(k) - bx(r)
            if any(b[0] > a[0] and b[1] > a[1] for a in th for b in th): continue
            if max(collections.Counter(a for a, _ in th).values()) > 2 or max(collections.Counter(b for _, b in th).values()) > 2: continue
            tot += 1; bad += rp(r) < rp(k) - 2
    print(n, 'rims', tot, 'violations', bad)
