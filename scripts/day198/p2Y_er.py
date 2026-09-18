"""
Day 198: Compute p_2(Y) . e_r(X) at m=4 for r=1,2,3 and expand in e-basis.

p_2(Y) = Y_1^2 + Y_2^2 + ... + Y_m^2 (power sum).

Framework copied from /home/agent/projects/proofs/scripts/day196/lp_test_length3.py.
"""
import sympy as sp
from itertools import combinations

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
    # Pad to equal length
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


def main():
    print("=" * 72)
    print("Day 198: p_2(Y) . e_r(X) for r=1,2,3")
    print("=" * 72)

    for r in [1, 2, 3]:
        n = r + 2  # degree of result
        m = max(4, n)  # need m >= n for the e-basis solve
        X, Ti_apply, Ti_inv_apply, Pi_apply, Y_apply = build_action(m)

        print()
        print("-" * 72)
        print(f"r = {r}  (degree of p_2(Y).e_r is n = {n}, using m = {m})")
        print("-" * 72)

        erX = e_r_X(m, r)
        # Compute p_2(Y) . e_r = sum_i Y_i^2 . e_r
        total = sp.Integer(0)
        for i in range(1, m+1):
            Yi_er = Y_apply(erX, i)
            Yi2_er = Y_apply(Yi_er, i)
            total = sp.expand(total + Yi2_er)

        print(f"p_2(Y) . e_{r}: #terms = {len(sp.Add.make_args(total))}")

        # Expand in e-basis of degree n
        expansion = expand_symmetric_in_e_basis(total, m, n)

        print()
        print(f"e-basis expansion of p_2(Y) . e_{r} (degree {n}):")

        # Sort partitions in reverse dominance-friendly order: by length ascending, then lex desc
        parts = partitions_of(n)
        for lam in parts:
            c = expansion.get(lam, sp.Integer(0))
            cs = sp.factor(sp.simplify(c))
            length = len(lam)
            print(f"  e_{lam}  (ell={length}): {cs}")

        # DS-triangular check: is support in {mu : mu >= (r,1,1) in dominance}?
        # For r=1: (r,1,1) = (1,1,1), partitions of 3. Only (3),(2,1),(1,1,1) all >= (1,1,1)?
        # Actually dominance of (1,1,1) means everyone >= (1,1,1) is everyone (since (1,1,1) is min).
        # More interesting for r=2: (r,1,1)=(2,1,1), partitions of 4. (4),(3,1),(2,2),(2,1,1) vs (1,1,1,1).
        # For r=3: (r,1,1)=(3,1,1), partitions of 5. Support should exclude (2,2,1),(2,1,1,1),(1^5)?
        # Let me just report which support is nonzero and check the (r,1,1) bound.
        print()
        pivot = tuple([r] + [1, 1])
        print(f"  DS-triangularity pivot (r,1,1) = {pivot}")
        print(f"  Support (nonzero coeffs) vs dominance >= {pivot}:")
        for lam in parts:
            c = expansion.get(lam, sp.Integer(0))
            cs = sp.simplify(c)
            is_zero = (cs == 0)
            ge_pivot = dominance_ge(lam, pivot)
            marker = ""
            if not is_zero and not ge_pivot:
                marker = "  <-- VIOLATES DS-triangularity (nonzero but NOT >= pivot)"
            if is_zero and ge_pivot:
                marker = "  (zero but >= pivot -- fine, DS is only upper bound)"
            print(f"    e_{lam}: zero={is_zero}, dominance>={pivot}? {ge_pivot}{marker}")


if __name__ == "__main__":
    main()
