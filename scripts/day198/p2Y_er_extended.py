"""
Day 198 extended: Compute p_2(Y) . e_r(X) at r=4 (m=6) and r=5 (m=7).

For each r, expand p_2(Y).e_r(X) in the e-basis of Lambda^(m) at degree n=r+2,
and check specific coefficients:
  - coef of e_(r+1,1) =?= (q-1)(q+1)/q^3 = (q^2-1)/q^3
  - coef of e_(r,2)   =?= -(q*t - q + t + 1)/q^3
  - coef of e_(r,1,1) =?= 1/q^3
  - all coeffs on partitions NOT dominating (r,1,1) are zero.

Framework copied from p2Y_er.py.
"""
import sympy as sp
from itertools import combinations
import time

q, t = sp.symbols('q t')


def qint(n):
    return sum(t**i for i in range(n))


def build_action(m):
    X = sp.symbols(f'X1:{m+1}')

    def si_apply(F, i):
        if i < 1 or i >= m:
            raise ValueError(f"si_apply supports i=1..m-1, got i={i}")
        subs = {X[i-1]: X[i], X[i]: X[i-1]}
        return sp.expand(F.xreplace(subs))

    def Ti_apply(F, i):
        F = sp.expand(F)
        sF = si_apply(F, i)
        diff = sp.expand(F - sF)
        quot = sp.cancel(diff / (X[i-1] - X[i]))
        result = t * sF + (t - 1) * (-X[i]) * quot
        return sp.expand(result)

    def Ti_inv_apply(F, i):
        Ti_F = Ti_apply(F, i)
        return sp.expand(Ti_F / t - (t - 1) / t * F)

    def Pi_apply(F):
        subs = {X[i]: X[i+1] for i in range(m-1)}
        subs[X[m-1]] = q**(-1) * X[0]
        Fshift = sp.expand(F.xreplace(subs))
        return sp.expand(X[0] * Fshift)

    def Y_apply(F, i):
        G = F
        for j in range(i, m):
            G = Ti_inv_apply(G, j)
        G = Pi_apply(G)
        for j in range(1, i):
            G = Ti_apply(G, j)
        return sp.expand(t**(m - i) * G)

    return X, Ti_apply, Ti_inv_apply, Pi_apply, Y_apply


def e_r_X(m, r):
    X = sp.symbols(f'X1:{m+1}')
    if r == 0:
        return sp.Integer(1)
    if r > m:
        return sp.Integer(0)
    result = sp.Integer(0)
    for combo in combinations(X, r):
        term = sp.Integer(1)
        for v in combo:
            term *= v
        result += term
    return sp.expand(result)


def partitions_of(n):
    result = []
    def rec(remaining, max_part, current):
        if remaining == 0:
            result.append(tuple(current))
            return
        for p in range(min(remaining, max_part), 0, -1):
            current.append(p)
            rec(remaining - p, p, current)
            current.pop()
    rec(n, n, [])
    return result


def e_lambda_X(m, lam):
    result = sp.Integer(1)
    for p in lam:
        result *= e_r_X(m, p)
    return sp.expand(result)


def expand_symmetric_in_e_basis(F, m, n):
    X = sp.symbols(f'X1:{m+1}')
    parts = partitions_of(n)
    if m < n:
        raise ValueError("Need m >= n.")
    canonicals = []
    for lam in parts:
        mon = sp.Integer(1)
        for i, p in enumerate(lam):
            mon *= X[i]**p
        canonicals.append(mon)
    Fpoly = sp.Poly(F, *X)
    b_vec = []
    for mon in canonicals:
        b_vec.append(Fpoly.coeff_monomial(sp.Poly(mon, *X).monoms()[0]))
    b_vec = sp.Matrix(b_vec)
    A_rows = []
    for lam in parts:
        el = e_lambda_X(m, lam)
        el_poly = sp.Poly(el, *X)
        row = [el_poly.coeff_monomial(sp.Poly(mon, *X).monoms()[0]) for mon in canonicals]
        A_rows.append(row)
    A = sp.Matrix(A_rows).T
    coeffs = A.solve(b_vec)
    result = {}
    for i, lam in enumerate(parts):
        c = sp.cancel(coeffs[i])
        result[lam] = c
    return result


def dominance_ge(mu, nu):
    """Return True if mu >= nu in dominance order (both partitions of same n)."""
    L = max(len(mu), len(nu))
    mu_p = list(mu) + [0] * (L - len(mu))
    nu_p = list(nu) + [0] * (L - len(nu))
    s_mu = 0
    s_nu = 0
    for k in range(L):
        s_mu += mu_p[k]
        s_nu += nu_p[k]
        if s_mu < s_nu:
            return False
    return True


def run_case(r, m):
    n = r + 2
    print()
    print("=" * 72)
    print(f"r = {r},  m = {m},  degree n = r+2 = {n}")
    print("=" * 72)
    t0 = time.time()

    X, Ti_apply, Ti_inv_apply, Pi_apply, Y_apply = build_action(m)

    print(f"[t={time.time()-t0:.1f}s] Building e_{r}(X) in {m} variables...")
    erX = e_r_X(m, r)

    # Compute p_2(Y) . e_r = sum_i Y_i . (Y_i . e_r)
    print(f"[t={time.time()-t0:.1f}s] Computing p_2(Y) . e_{r} = sum_i Y_i^2 . e_{r}...")
    total = sp.Integer(0)
    for i in range(1, m+1):
        ti = time.time()
        Yi_er = Y_apply(erX, i)
        Yi2_er = Y_apply(Yi_er, i)
        total = sp.expand(total + Yi2_er)
        print(f"    [t={time.time()-t0:.1f}s]  Y_{i}^2 . e_{r}: cumulative #terms = {len(sp.Add.make_args(total))}  ({time.time()-ti:.1f}s for this i)")

    print(f"[t={time.time()-t0:.1f}s] Total p_2(Y).e_{r}: #terms = {len(sp.Add.make_args(total))}")

    # Expand in e-basis
    print(f"[t={time.time()-t0:.1f}s] Expanding in e-basis of degree {n}...")
    expansion = expand_symmetric_in_e_basis(total, m, n)

    print(f"[t={time.time()-t0:.1f}s] Factoring coefficients...")
    parts = partitions_of(n)
    print()
    print(f"FULL e-basis expansion of p_2(Y) . e_{r}  (degree {n}):")
    print("-" * 72)
    factored = {}
    for lam in parts:
        c = expansion.get(lam, sp.Integer(0))
        cs = sp.factor(sp.simplify(c))
        factored[lam] = cs
        length = len(lam)
        print(f"  e_{lam}  (ell={length}): {cs}")

    # Specific checks
    print()
    print("SPECIFIC COEFFICIENT CHECKS:")
    print("-" * 72)

    pivot = (r, 1, 1)

    # Check 1: coef of e_(r+1,1)
    lam1 = (r+1, 1)
    expected1 = (q - 1) * (q + 1) / q**3
    got1 = factored.get(lam1, sp.Integer(0))
    diff1 = sp.simplify(sp.together(got1 - expected1))
    print(f"  coef of e_{lam1}: {got1}")
    print(f"    expected (q-1)(q+1)/q^3 = {sp.factor(expected1)}")
    print(f"    match? {diff1 == 0}")

    # Check 2: coef of e_(r,2)
    lam2 = (r, 2)
    expected2 = -(q*t - q + t + 1) / q**3
    got2 = factored.get(lam2, sp.Integer(0))
    diff2 = sp.simplify(sp.together(got2 - expected2))
    print(f"  coef of e_{lam2}: {got2}")
    print(f"    expected -(q*t - q + t + 1)/q^3 = {sp.factor(expected2)}")
    print(f"    match? {diff2 == 0}")

    # Check 3: coef of e_(r,1,1)
    lam3 = (r, 1, 1)
    expected3 = sp.Integer(1) / q**3
    got3 = factored.get(lam3, sp.Integer(0))
    diff3 = sp.simplify(sp.together(got3 - expected3))
    print(f"  coef of e_{lam3}: {got3}")
    print(f"    expected 1/q^3 = {sp.factor(expected3)}")
    print(f"    match? {diff3 == 0}")

    # Check 4: all coeffs on partitions NOT dominating (r,1,1) are zero.
    print()
    print(f"  Support restriction check: all coefs on lambda NOT >= (r,1,1) should be zero.")
    all_ok = True
    for lam in parts:
        c = factored.get(lam, sp.Integer(0))
        cs = sp.simplify(c)
        is_zero = (cs == 0)
        ge_pivot = dominance_ge(lam, pivot)
        if not ge_pivot:
            # This partition is NOT >= (r,1,1); its coefficient must be zero.
            status = "OK (zero)" if is_zero else f"VIOLATION (nonzero: {cs})"
            if not is_zero:
                all_ok = False
            print(f"    e_{lam} (NOT >= pivot): {status}")
    print(f"  All non-dominating coeffs zero? {all_ok}")

    print()
    print(f"[t={time.time()-t0:.1f}s] r={r} complete.")
    return factored


def main():
    print("=" * 72)
    print("Day 198 extended: p_2(Y) . e_r(X) for r=4, r=5")
    print("=" * 72)

    # r = 4, m = 6
    factored_r4 = run_case(r=4, m=6)

    # r = 5, m = 7 (only if r=4 didn't blow up)
    factored_r5 = run_case(r=5, m=7)


if __name__ == "__main__":
    main()
