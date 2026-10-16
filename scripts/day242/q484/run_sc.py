"""Structure constants, support bound, Cor 10(2), Laurent census vs Clio's jump criterion."""
import sys, sympy as sp, json
from engine import *
X, N, qarg, norm = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4] == 'norm'
qq = q if qarg == 'q' else sp.Rational(qarg)
A, C, w = structure_constants(X, N, qq, normalise=norm)
S = subsets(N)
sup_ok = sup_tot = 0; cor_ok = cor_tot = 0
census = {'laurent': 0, 'pole': 0, 'agree': 0, 'FP': 0, 'FN': 0, 'tot': 0}
witness = None
for (l, m), c in C.items():
    for nu, val in c.items():
        allowed = geq(nu, l) and geq(nu, m)
        if not allowed:
            sup_tot += 1; sup_ok += (val == 0)
            continue
        census['tot'] += 1
        L = is_laurent(val, w)
        census['laurent' if L else 'pole'] += 1
        pred = (len(nu) - max(len(l), len(m))) <= 1   # Clio Conj 15 uses |.| = cardinality? see below
        predw = (sum(nu) - max(sum(l), sum(m))) <= 1
        census.setdefault('agree_card', 0); census.setdefault('agree_wt', 0)
        census['agree_card'] += (pred == L); census['agree_wt'] += (predw == L)
        if l == () and m == () and nu == tuple(range(1, N+1)) and N == 2:
            witness = sp.factor(val)
    if geq(l, m):
        cor_tot += 1; cor_ok += sp.cancel(c[l] - A[l, m]) == 0
print(X, N, qarg, 'norm' if norm else 'unnorm', 'support-zero', sup_ok, '/', sup_tot,
      'Cor10(2)', cor_ok, '/', cor_tot, census, flush=True)
if witness is not None: print('  c^{(1,2)}_{00} =', witness, flush=True)
