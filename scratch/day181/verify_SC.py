"""
Day 181 — Verify sub-claim (SC) at n=5, 6 with focus on m'' containing E_3.

(SC): For r >= 1 and m'' in Q[E_1, E_2, E_3]:
      T^{X,r}(m'') mod E_{>=4} has rho <= 2r + rho(m'').

Where:
    T^{X,r}(m'') = sum_{i<j} (u_i + u_j + 1) * X_{ij}^r * m''|_{ij}
    X_{ij}       = E_3|_{ij} - E_3
                 = (2 E_2 + E_1) - (u_i + u_j)(E_1 + 1) + (u_i^2 + u_j^2)
    m''|_{ij}    = m''(u + e_i + e_j)  (shift u_i, u_j by 1).

We first re-verify §4.3 for (r=1, m''=1) and (r=1, m''=E_2), then check
cases involving E_3 in m''.
"""

import sympy as sp
from itertools import combinations

def rho(k):
    return (k + 1) // 2

def do_verification(n):
    print(f"===== n = {n} =====")
    u = sp.symbols(f'u1:{n+1}')
    E_syms = sp.symbols(f'E1:{n+1}')

    def E_polys():
        return [sp.symmetric_poly(k, u) for k in range(1, n+1)]

    E_vals = E_polys()

    def shift_ij(expr, i, j):
        return expr.subs([(u[i-1], u[i-1]+1), (u[j-1], u[j-1]+1)], simultaneous=True)

    def X_ij(i, j):
        # X_{ij} = E_3|_{ij} - E_3
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
            return -1  # sentinel for zero
        p = sp.Poly(sp.expand(expr_in_E), *E_syms)
        return max(rho_weight_of_monom(monom) for monom, _ in p.terms())

    def m_from_indices(indices):
        prod = sp.Integer(1)
        for k in indices:
            prod *= E_vals[k-1]
        return sp.expand(prod)

    def rho_of_indices(indices):
        return sum(rho(k) for k in indices)

    # Test cases: (name, indices for m'', r)
    tests = []
    # r = 1 cases
    for r in [1, 2]:
        tests += [
            ("1", [], r),
            ("E_1", [1], r),
            ("E_2", [2], r),
            ("E_3", [3], r),
            ("E_1*E_3", [1, 3], r),
            ("E_2*E_3", [2, 3], r),
            ("E_3^2", [3, 3], r),
            ("E_1^2*E_3", [1, 1, 3], r),
            ("E_2^2", [2, 2], r),
        ]

    records = []
    for name, indices, r in tests:
        rho_m = rho_of_indices(indices)
        target_bound = 2*r + rho_m  # (SC) says rho <= this
        m_expr = m_from_indices(indices)
        try:
            TXr = T_Xr(m_expr, r)
            TXr_E = symmetric_to_E(TXr)
            TXr_reduced = drop_E_geq_4(TXr_E)
            mr = max_rho(TXr_reduced)
            ok = (mr <= target_bound)
            zero_flag = "(zero)" if mr == -1 else ""
            print(f"  r={r} m''={name:12s} rho(m'')={rho_m}  bound={target_bound}  observed max_rho={mr:3d}  {'PASS' if ok else 'FAIL'} {zero_flag}")
            records.append((name, r, rho_m, target_bound, mr, ok))
        except Exception as e:
            print(f"  ERROR for r={r} m''={name}: {e}")

    passed = sum(1 for r in records if r[5])
    print(f"  TOTAL: {passed}/{len(records)} pass")
    return records

if __name__ == "__main__":
    do_verification(5)
    do_verification(6)
