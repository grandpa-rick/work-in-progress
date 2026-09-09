"""
Day 181 — Verify the u_i u_j-lemma: sum_{i<j} (u_i^a u_j^b + u_i^b u_j^a)
has top-rho <= (a+b) - 1 mod E_{>=4}, for a, b >= 1.

More precisely: for any symmetric polynomial f(u_i, u_j) divisible by u_i u_j
(equivalently, f(u_i, u_j) = u_i u_j g(u_i, u_j) with g symmetric), the sum
sum_{i<j} f has top rho <= u-deg(f) - 1 mod E_{>=4}.
"""

import sympy as sp
from itertools import combinations

def rho(k):
    return (k + 1) // 2

def do_check(n=5):
    print(f"n = {n}")
    u = sp.symbols(f'u1:{n+1}')
    E_syms = sp.symbols(f'E1:{n+1}')

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

    tests = []
    # (a, b) with a, b >= 1
    for a in range(1, 4):
        for b in range(1, 4):
            tests.append((a, b))

    all_pass = True
    for (a, b) in tests:
        # sum_{i<j} (u_i^a u_j^b + u_i^b u_j^a)
        S = sum(u[i]**a * u[j]**b + u[i]**b * u[j]**a for i, j in combinations(range(n), 2))
        S_E = symmetric_to_E(S)
        target = a + b  # naive top rho
        top = top_rho_slice(S_E, target)
        ok = sp.expand(top) == 0
        print(f"  a={a} b={b}: sum drops top rho={target}?  top slice = {top}  {'PASS' if ok else 'FAIL'}")
        if not ok:
            all_pass = False
    print()
    print(f"Overall: {'ALL PASS' if all_pass else 'SOME FAIL'}")
    return all_pass

if __name__ == "__main__":
    do_check(n=5)
