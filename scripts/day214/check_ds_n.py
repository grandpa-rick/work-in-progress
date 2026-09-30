# parallel per-partition version of check_ds_all_lengths for a single size n
import sympy as sp, sys
from multiprocessing import Pool
from ek_subset_engine import *
S, T = sp.Rational(3,7), sp.Rational(-5,11)
def job(lam):
    n = sum(lam); xs = setup(n)
    out = e_expand(e_lam_star(lam, xs, S, T), xs)
    supp = all(dominates(mu, lam) for mu in out); lead = out.get(lam, 0) == S**nstat(lam)
    full = set(out) == {mu for mu in parts(n) if dominates(mu, lam)}
    V = e_expand(e_lam_star(lam, xs, sp.Integer(1), T), xs) == {lam: 1}
    return lam, supp, lead, V, full
if __name__ == '__main__':
    n = int(sys.argv[1]); ok = True
    with Pool(6) as P:
        for lam, supp, lead, V, full in P.imap(job, list(parts(n))):
            ok &= supp and lead and V
            print(lam, 'supp', supp, 'lead', lead, 's=1', V, 'full-upset', full, flush=True)
    print('ALL OK' if ok else 'FAIL')
