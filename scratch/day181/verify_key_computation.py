"""
Day 181 — Verify the KEY COMPUTATION in the (SC) proof.

The proof rests on: M_l(f)^{[top,rho=2r+l]} = 0 for r >= 1,
where f(u) = (A - B*u + u^2)^r, A = 2*E_2 + E_1, B = E_1 + 1.

Since M_l(f) = sum_i u_i^l f(u_i) = sum_{b+d=r} binom(r,b) (-1)^b B^b p_{l+b+2d},
top rho piece = E_1^l * E_1^{2r} * (1-1)^r = 0 for r >= 1.

Verify this via direct sympy computation for r = 2, 3, 4 at n = 5.
"""

import sympy as sp
from itertools import combinations

def rho(k):
    return (k + 1) // 2

def do_check(n=5, l_test=[0, 1, 2], r_test=[1, 2, 3, 4]):
    print(f"n = {n}")
    u = sp.symbols(f'u1:{n+1}')
    E_syms = sp.symbols(f'E1:{n+1}')
    E_vals = [sp.symmetric_poly(k, u) for k in range(1, n+1)]

    A_val = 2 * E_vals[1] + E_vals[0]  # 2*E_2 + E_1
    B_val = E_vals[0] + 1  # E_1 + 1

    def M_l(f_expr, l):
        return sum(u[i]**l * f_expr.subs(sp.Symbol('u_dummy'), u[i]) for i in range(n))

    def symmetric_to_E(expr):
        p_c = sp.expand(expr)
        result, rem, _ = sp.symmetrize(p_c, u, formal=True, symbols=E_syms)
        if sp.expand(rem) != 0:
            raise ValueError(f"symmetrize left remainder: {rem}")
        return sp.expand(result)

    def rho_weight_of_monom(monom):
        return sum(a * rho(k) for k, a in enumerate(monom, start=1))

    def top_rho_slice(expr_in_E, target_rho):
        p = sp.Poly(sp.expand(expr_in_E), *E_syms)
        result = sp.Integer(0)
        for monom, coeff in p.terms():
            # Drop E_{>=4}
            if any(monom[k-1] > 0 for k in range(4, n+1)):
                continue
            w = rho_weight_of_monom(monom)
            if w != target_rho:
                continue
            term = coeff
            for k, a in enumerate(monom, start=1):
                term *= E_syms[k-1]**a
            result += term
        return sp.expand(result)

    def f_of_u(r):
        udum = sp.Symbol('u_dummy')
        return (A_val - B_val * udum + udum**2)**r

    all_pass = True
    for r in r_test:
        fu = f_of_u(r)
        for l in l_test:
            udum = sp.Symbol('u_dummy')
            M = sum(u[i]**l * fu.subs(udum, u[i]) for i in range(n))
            M_E = symmetric_to_E(M)
            top = top_rho_slice(M_E, 2*r + l)
            ok = sp.expand(top) == 0
            print(f"  r={r} l={l}: target rho={2*r+l}  M_l(f)^[top] = {top}  {'PASS (0)' if ok else 'FAIL'}")
            if not ok:
                all_pass = False
    print()
    print(f"Overall: {'ALL PASS' if all_pass else 'SOME FAIL'}")
    return all_pass

if __name__ == "__main__":
    do_check(n=5)
