"""Day 210: extract e-basis coefficients of e_k*(e_a e_b) = t^{-C(k,2)} e_k(Y).(e_a e_b)
directly from the statement of (TC) (Day 209 §0), written independently of Day 208 code.
T_k = sum_{b0,i,j} s^{2b0} t^{-b0(i+j)} V^{(k-b0)}_{ij} e_{b0} E(t^i z) E(t^j w)."""
import sympy as sp
from functools import lru_cache
s, t, z, w = sp.symbols('s t z w')

def qp(a, n):  # (a;t)_n
    r = sp.Integer(1)
    for i in range(n): r *= (1 - a*t**i)
    return r
@lru_cache(None)
def alpha(j):
    if j < 0: return sp.Integer(0)
    r = sp.Integer(1)
    for i in range(1, j+1): r *= (s - t**i)
    return r/qp(t, j)
@lru_cache(None)
def c(n, j):
    if j < 0 or j > n: return sp.Integer(0)
    return qp(s, n-j)/qp(t, n-j)*(alpha(j) - s*t**(n-j)*alpha(j-1))
def N(n, j): return t**(-n*j)*c(n, j)
def K(i, j):
    r = sp.Integer(1)
    for p in range(i): r *= (t**p*z - s*w)/(t**p*z - t**j*w)
    for q in range(j): r *= (s*z - t**q*w)/(t**i*z - t**q*w)
    return r

def coeffs(k, a, b):
    """dict partition -> coefficient (in s,t) of e_mu in e_k*(e_a e_b), unique expansion in Lambda."""
    d = a + b
    acc = {}
    for b0 in range(k+1):
        n = k - b0
        for i in range(n+1):
            for j in range(n+1-i):
                Kij = K(i, j)
                for n1 in range(i, n-j+1):
                    n2 = n - n1
                    pref = (s**(2*b0)*t**(-b0*(i+j)) * s**((n1-i)+(n2-j))
                            * t**(-(n1-i)*j-(n2-j)*i) * N(n1, i)*N(n2, j))
                    if pref == 0: continue
                    tot = d + n1 + n2
                    for r1 in range(tot+1):
                        r2 = tot - r1
                        mu = tuple(sorted([x for x in (b0, r1, r2) if x > 0], reverse=True))
                        term = pref*Kij*t**(i*r1 + j*r2)*z**(r1-n1)*w**(r2-n2)
                        acc[mu] = acc.get(mu, 0) + term
    out = {}
    for mu, R in acc.items():
        R = sp.factor(sp.together(R.subs(w, 1)))
        num, den = sp.fraction(R)
        assert sp.Poly(den, z).degree() == 0, (mu, den)   # poles cancel (TC Corollary)
        cf = sp.factor(sp.Poly(sp.expand(num), z).coeff_monomial(z**a)/den)
        if cf != 0: out[mu] = cf
    return out

def ek_er(k, r):
    """207b: e_k*e_r = sum_b s^b F_{k-b}(t^{r-b}) e_b e_{r+k-b}."""
    out = {}
    for b0 in range(k+1):
        n = k - b0
        F = sum(c(n, j)*t**((r-b0)*j) for j in range(n+1))
        mu = tuple(sorted([x for x in (b0, r+k-b0) if x > 0], reverse=True))
        out[mu] = out.get(mu, 0) + s**b0*F
    return {m: sp.factor(v) for m, v in out.items() if sp.simplify(v) != 0}
