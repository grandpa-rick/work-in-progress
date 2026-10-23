# Fit b_k(j)=H_k(j)-H_k(j-1) = sum_{rho<k} c_rho H_rho(j) over rho = k minus <=2 boxes per removable row-end, report uniqueness.
import sys, itertools, sympy as sp
N_ARG = sys.argv[1] if len(sys.argv) > 1 else '5'; DEP = sys.argv[2:]; sys.argv = sys.argv[:1] + [N_ARG]
import gate as G
from fractions import Fraction as Fr
x = sp.Symbol('j')
def poly(k):
    v = G.Hvals(k); return sp.expand(sp.interpolate(list(zip(range(1 - len(v), 1), v)), x)) if len(v) > 1 else sp.Integer(v[0])
def fit(k, depth):
    sup = [r for r in G.below(k) if r != k and sum(k) - sum(r) <= depth]
    cs = sp.symbols(f'c0:{len(sup)}'); b = sp.expand(poly(k) - poly(k).subs(x, x - 1))
    e = sp.Poly(b - sum(c * poly(r) for c, r in zip(cs, sup)), x).all_coeffs()
    sol = sp.solve(e, cs, dict=True)
    return sup, cs, sol
if __name__ == '__main__':
    for k in [tuple(sorted(c, reverse=True)) for L in range(1, G.N + 1) for c in itertools.combinations(range(1, G.N + 1), L)]:
        sup, cs, sol = fit(k, int(DEP[0]) if DEP else 4)
        print(k, {''.join(map(str, r)) or '0': sol[0].get(c, c) for c, r in zip(cs, sup)} if sol else None)
