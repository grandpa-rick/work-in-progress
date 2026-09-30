import sys; sys.path.insert(0,'../../proofs/scripts/day196')
import sympy as sp, time
from ds_test_221 import compute_e_star_product, expand_symmetric_in_e_basis
from ds3 import elam
q,t = sp.symbols('q t'); s = sp.Symbol('s')
for lam, m in [((2,2,1),5), ((3,2,1),6)]:
    t0=time.time()
    F = compute_e_star_product(m, list(lam))
    ex = expand_symmetric_in_e_basis(F, m, sum(lam))
    tc = elam((lam[2], lam[0], lam[1]))
    bad = [mu for mu in set(ex)|set(tc) if sp.cancel(sp.sympify(ex.get(mu,0)) - sp.sympify(tc.get(mu,0)).subs(s,1/q)) != 0]
    print(lam, 'm=',m, 'direct AHA vs TC-route:', 'MATCH' if not bad else f'MISMATCH {bad}', f'{time.time()-t0:.0f}s', flush=True)
