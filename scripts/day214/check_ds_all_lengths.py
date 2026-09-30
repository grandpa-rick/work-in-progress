import sympy as sp, sys, pickle, random
from ek_subset_engine import *
S, T = sp.Rational(3,7), sp.Rational(-5,11)
ok = True
# 1) cross-check vs Day 211 (TC-based) engine at numeric point
st = pickle.load(open('../day211/ds3_full_N6.pkl','rb')) if __import__('os').path.exists('../day211/ds3_full_N6.pkl') else None
NMAX = int(sys.argv[1])
for n in range(1, NMAX+1):
    xs = setup(n)
    for lam in parts(n):
        out = e_expand(e_lam_star(lam, xs, S, T), xs)
        supp = all(dominates(mu, lam) for mu in out)
        lead = out.get(lam, 0) == S**nstat(lam)
        up = [mu for mu in parts(n) if dominates(mu, lam)]
        full = set(out) == set(up)
        out1 = e_expand(e_lam_star(lam, xs, sp.Integer(1), T), xs)
        V = (out1 == {lam: 1})
        # order independence (commutativity) sanity
        out2 = e_expand(e_lam_star(lam, xs, S, T, order=list(lam)), xs)
        comm = (out2 == out)
        cross = ''
        if st is not None and len(lam) == 3 and lam in st:
            s_, t_ = sp.symbols('s t')
            ref = {mu: sp.nsimplify(sp.cancel(c.subs({sp.Symbol('s'): S, sp.Symbol('t'): T}))) for mu, c in st[lam].items()}
            ref = {mu: c for mu, c in ref.items() if c != 0}
            cross = f' vsTC:{ref == out}'
            ok &= ref == out
        ok &= supp and lead and V and comm
        print(lam, 'supp', supp, 'lead', lead, 's=1', V, 'comm', comm, 'full-upset', full, cross, flush=True)
print('ALL OK' if ok else 'FAIL')
