"""
Extend (SC) verification to n=7 for a subset of cases (r=1 only, small m'').
"""

import sympy as sp
from itertools import combinations

def rho(k):
    return (k + 1) // 2

def do_verification(n):
    print(f"===== n = {n} =====")
    u = sp.symbols(f'u1:{n+1}')
    E_syms = sp.symbols(f'E1:{n+1}')
    E_vals = [sp.symmetric_poly(k, u) for k in range(1, n+1)]

    def shift_ij(expr, i, j):
        return expr.subs([(u[i-1], u[i-1]+1), (u[j-1], u[j-1]+1)], simultaneous=True)

    def X_ij(i, j):
        return sp.expand(shift_ij(E_vals[2], i, j) - E_vals[2])

    def T_Xr(m_expr, r):
        idx = list(range(1, n+1))
        total = sp.Integer(0)
        for i, j in combinations(idx, 2):
            mij = shift_ij(m_expr, i, j)
            weight = u[i-1] + u[j-1] + 1
            Xij = X_ij(i, j)
            total += weight * (Xij ** r) * mij
        return sp.expand(total)

    def symmetric_to_E(poly_expr):
        p_c = sp.expand(poly_expr)
        result, rem, _ = sp.symmetrize(p_c, u, formal=True, symbols=E_syms)
        if sp.expand(rem) != 0:
            raise ValueError(f"symmetrize left remainder: {rem}")
        return sp.expand(result)

    def rho_weight_of_monom(monom):
        return sum(a * rho(k) for k, a in enumerate(monom, start=1))

    def drop_E_geq_4(expr_in_E):
        p = sp.Poly(sp.expand(expr_in_E), *E_syms)
        result = sp.Integer(0)
        for monom, coeff in p.terms():
            if any(monom[k-1] > 0 for k in range(4, n+1)):
                continue
            term = coeff
            for k, a in enumerate(monom, start=1):
                term *= E_syms[k-1]**a
            result += term
        return sp.expand(result)

    def max_rho(expr_in_E):
        if sp.expand(expr_in_E) == 0:
            return -1
        p = sp.Poly(sp.expand(expr_in_E), *E_syms)
        return max(rho_weight_of_monom(monom) for monom, _ in p.terms())

    def m_from_indices(indices):
        prod = sp.Integer(1)
        for k in indices:
            prod *= E_vals[k-1]
        return sp.expand(prod)

    def rho_of_indices(indices):
        return sum(rho(k) for k in indices)

    # Reduced test set for n=7 (compute is slow)
    tests = [
        ("1", [], 1),
        ("E_1", [1], 1),
        ("E_3", [3], 1),
        ("E_2", [2], 1),
        ("E_1*E_3", [1, 3], 1),
    ]

    records = []
    for name, indices, r in tests:
        rho_m = rho_of_indices(indices)
        target_bound = 2*r + rho_m
        m_expr = m_from_indices(indices)
        try:
            TXr = T_Xr(m_expr, r)
            TXr_E = symmetric_to_E(TXr)
            TXr_reduced = drop_E_geq_4(TXr_E)
            mr = max_rho(TXr_reduced)
            ok = (mr <= target_bound)
            print(f"  r={r} m''={name:12s} rho(m'')={rho_m}  bound={target_bound}  observed max_rho={mr:3d}  {'PASS' if ok else 'FAIL'}")
            records.append((name, r, rho_m, target_bound, mr, ok))
        except Exception as e:
            print(f"  ERROR for r={r} m''={name}: {e}")

    passed = sum(1 for rrr in records if rrr[5])
    print(f"  TOTAL: {passed}/{len(records)} pass")
    return records

if __name__ == "__main__":
    do_verification(7)
