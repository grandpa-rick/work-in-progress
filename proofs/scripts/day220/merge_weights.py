"""Day 220: merge weights W_k(J) = lin_e [[E_k^{(p)}, e_j1],...,e_jp](1), p=|J|, n=k+|J|.
lin_e G = (-1)^{n-1} G(1,z,...,z^{n-1}), z primitive n-th root (all e_mu, mu != (n), vanish there).
Top symbol: [[E_k^{(p)},a_1..a_p]](1) = sum_A c_A x_A prod Delta_A(a_i); at zeta-point Delta_A e_j = (-1)^{j-1} p_j(X_A).
So W_k(J) = (-1)^{n-1+sum(j-1)} sum_{|A|=k} c_A x_A prod_i p_{j_i}(X_A) at x=zeta-point.  Grade: computed."""
import sys, itertools
from sympy import symbols, Poly, cyclotomic_poly, invert, factor, expand, rem, QQ
z, t = symbols('z t')
def W(k, J):
    n = k + sum(J); Phi = cyclotomic_poly(n, z)
    tot = 0
    for A in itertools.combinations(range(n), k):
        B = [j for j in range(n) if j not in A]
        num = 1; den = 1
        for i in A:
            for j in B:
                num = expand(num*(z**i - t*z**j)); den = rem(expand(den*(z**i - z**j)), Phi, z)
        term = num*z**sum(A)
        for jj in J: term = expand(term*sum(z**(jj*i) for i in A))
        term = rem(expand(term*invert(den, Phi, z)), Phi, z)
        tot = rem(expand(tot + term), Phi, z)
    tot = Poly(tot, z)
    assert tot.degree() <= 0, tot
    val = tot.as_expr()
    sgn = (-1)**(n-1+sum(j-1 for j in J))
    return factor(sgn*val)
def parts(n, m=None):
    if m is None: m = n
    if n == 0: yield (); return
    for k in range(min(n, m), 0, -1):
        for p in parts(n-k, k): yield (k,)+p
if __name__ == '__main__':
    MAXN = int(sys.argv[1])
    for n in range(2, MAXN+1):
        for k in range(1, n):
            for J in parts(n-k):
                w = W(k, J); p = len(J)
                ok = all(c > 0 for c in Poly(expand((-1)**p*w), t).coeffs())
                print(f'n={n} W_{k}{J} = {w}   (-1)^p W in N[t]: {ok}', flush=True)
