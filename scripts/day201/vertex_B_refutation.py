"""
Day 201: Refute Vertex B (Jack shape check for p_2(Y) vs Thibon Delta_2(alpha)).

STRATEGY:
Rick's p_2(Y) is degree +2 on X (verified: p_2(Y).e_r has deg r+2).
Thibon's Delta_2(alpha) is degree 0.  Direct comparison impossible.

FALLBACK COMPARISON (spectral / ⋆-eigenvalue):
Compute F := p_2(Y).1 as a symmetric polynomial in X, then compare
its Macdonald-spectrum eigenvalue in the Jack limit against Thibon's
Delta_2(alpha) eigenvalue on Jack polynomials.

Rick's F = p_2(Y).1 = (p_2 - (t-1)(q+1) e_2) / q  (verified at m=3 SymPy).

Substitute q = t^alpha, t = 1 + eps, spec_lambda_i = q^{lambda_i} t^{-(i-1)}
= (1+eps)^{alpha*lambda_i - (i-1)} = 1 + eps*(alpha*lambda_i - (i-1)) + O(eps^2).

F(spec_lambda) = p_2(spec)/q - (t-1)(q+1) e_2(spec)/q.

Expand in eps and compare to Thibon 2 C_1^(alpha)(lambda) = alpha*(sum lambda_i^2 - |lambda|) - 2 n(lambda).
"""
import sympy as sp
from itertools import combinations

q, t, alpha, eps = sp.symbols('q t alpha epsilon')

# Rick's F(X) = (p_2 - (t-1)(q+1) e_2)/q where the symmetric functions
# are evaluated on the Macdonald spectrum spec_lambda_i = q^{lambda_i} t^{-(i-1)}.

def spec_i(lam, i, m):
    """spec_lambda_i for i in 1..m (padding lambda with zeros)."""
    lam_i = lam[i-1] if i <= len(lam) else 0
    return q**lam_i * t**(-(i-1))

def p_k_spec(lam, k, m):
    return sum(spec_i(lam, i, m)**k for i in range(1, m+1))

def e_2_spec(lam, m):
    result = sp.Integer(0)
    for i in range(1, m+1):
        for j in range(i+1, m+1):
            result += spec_i(lam, i, m) * spec_i(lam, j, m)
    return result

def F_spec(lam, m):
    """F(spec_lambda) = (p_2(spec) - (t-1)(q+1) e_2(spec))/q"""
    return (p_k_spec(lam, 2, m) - (t-1)*(q+1)*e_2_spec(lam, m)) / q

def n_stat(lam):
    return sum((i)*lam[i] for i in range(len(lam)))

def content_squared_sum(lam):
    """2 C_1^(alpha)(lam) = alpha (sum lam_i^2 - |lam|) - 2 n(lam) 
    where n(lam) = sum (i-1) lam_i.  This is Thibon's Delta_2(alpha) eigenvalue.
    """
    L = sum(l**2 for l in lam) - sum(lam)
    return alpha * L - 2 * n_stat(lam)

def jack_expand(expr, order):
    """Substitute q = (1+eps)^alpha, t = 1+eps, expand to O(eps^order)."""
    expr2 = expr.subs(q, (1+eps)**alpha).subs(t, 1+eps)
    ser = sp.series(expr2, eps, 0, order+1).removeO()
    return sp.expand(ser)

def compare_lambdas():
    print("="*72)
    print("Vertex B refutation: F(spec_lambda) vs Thibon's Delta_2(alpha) eigvals")
    print("="*72)
    print()
    m = 5  # enough for small lambdas
    test_lambdas = [(), (1,), (2,), (1,1), (3,), (2,1), (1,1,1), (2,2), (3,1)]
    
    print(f"m = {m} (num variables)")
    print()
    for lam in test_lambdas:
        F_val = F_spec(lam, m)
        F_series = jack_expand(F_val, 2)
        thibon_val = content_squared_sum(lam)
        
        # Extract eps^0 and eps^1 coefficients
        F_series_poly = sp.Poly(sp.expand(F_series), eps)
        F0 = sp.simplify(F_series_poly.coeff_monomial(eps**0))
        F1 = sp.simplify(F_series_poly.coeff_monomial(eps**1))
        
        print(f"  lambda = {lam}:")
        print(f"    F(spec)     = {F0} + eps*({F1}) + O(eps^2)")
        print(f"    Thibon eigval  = 2 C_1^(alpha) = {sp.simplify(thibon_val)}")
        # For Vertex B rescue: check if F(spec) - m corresponds to Thibon
        # after some rescale. Difference = F1 - c*thibon_val for various c.
        
        # Try F1 = 2 C_1 + constant? Constant should be m-independent, but F1 depends on m.
        # For lambda = () (empty): 2 C_1 = 0. F1 = ? -2m(m-1) - alpha*m
        # Just print eigenvalue structure differences

    print()
    print("STRUCTURAL DIAGNOSTIC")
    print("-" * 40)
    print("F(spec_lambda) is a degree-2 polynomial in spec_lambda_i's.")
    print("Each spec_i = q^{lambda_i} t^{-(i-1)} ~ 1 + eps*(alpha*lambda_i - (i-1)).")
    print("So F(spec) at O(eps^1) is LINEAR in lambda_i (not quadratic).")
    print("Thibon's 2 C_1^(alpha)(lambda) = alpha*(sum lambda_i^2 - |lambda|) - 2 n(lambda)")
    print("is QUADRATIC in lambda_i.")
    print()
    print("=> NO SCALAR RESCALING makes F(spec)/eps match Thibon's eigenvalue.")
    print("=> Vertex B REFUTED at leading order in eps.")
    print()
    print("Alternative reading: does F(spec) at O(eps^2) become quadratic")
    print("and match Thibon?  Let's test by computing eps^2 coefficient.")
    print()
    for lam in test_lambdas:
        F_val = F_spec(lam, m)
        F_series = jack_expand(F_val, 2)
        F_series_poly = sp.Poly(sp.expand(F_series), eps)
        F2 = sp.simplify(F_series_poly.coeff_monomial(eps**2))
        thibon_val = content_squared_sum(lam)
        
        # If F(spec) - m - eps*(linear stuff) = eps^2 * (quadratic in lambda), 
        # maybe (F2 - constants) = c(alpha) * thibon_val for some c(alpha)?
        # Diagnostic: F2 - F2|_lambda=empty, is this proportional to thibon_val?
        empty_F2 = sp.Poly(sp.expand(jack_expand(F_spec((), m), 2)), eps).coeff_monomial(eps**2)
        F2_shifted = sp.simplify(F2 - empty_F2)
        print(f"  lambda = {lam}: F(spec)@eps^2 - F(empty)@eps^2 = {F2_shifted}")
        print(f"                Thibon eigval = {sp.simplify(thibon_val)}")

if __name__ == "__main__":
    compare_lambdas()
