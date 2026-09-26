"""Day 207b: the E_k induction step, checked on closed forms.
T_k^{(m)}(X;z) := sum_{j,b} s^b t^{-kj} c(k-b,j) e_b(X) z^{b-k} E(X;t^j z).
Claim (proved via Lemma 2 + (*)): [k] T_k^{(m)}(X) = sum_i X_i(1+sX_i z) prod_{j!=i} a_ij T_{k-1}^{(m-1)}(X^_i).
Exact rationals for X,s,t; z symbolic (Laurent polynomial identity)."""
import sympy as sp, random, itertools
def _prod(L):
    r = sp.Integer(1)
    for v in L: r = r*v
    return r
from fractions import Fraction as Fr
z = sp.Symbol('z')
def cF(n, j, sv, tv):
    if j < 0 or n - j < 0: return sp.Integer(0)
    po = lambda a, N: _prod([1 - a*tv**i for i in range(N)])
    al = lambda J: sp.Integer(0) if J < 0 else _prod([sv - tv**i for i in range(1, J+1)])/po(tv, J)
    return po(sv, n-j)/po(tv, n-j)*(al(j) - sv*tv**(n-j)*al(j-1))
def e(X, b):
    if b < 0 or b > len(X): return 0
    return sum(_prod(c) for c in itertools.combinations(X, b))
def E(X, w): return _prod([1 + x*w for x in X])
def Tk(k, X, sv, tv):
    return sp.expand(sum(sv**b*tv**(-k*j)*cF(k-b, j, sv, tv)*e(X, b)*z**(b-k)*E(X, tv**j*z)
                         for b in range(k+1) for j in range(k-b+1)))
random.seed(5); ok = True
for m in range(1, 6):
    for k in range(1, 5):
        X = []
        while len(X) < m:
            v = sp.Rational(random.randint(1, 30), random.randint(31, 70))
            if v not in X: X.append(v)  # distinct (a_ij has X_i - X_j in the denominator)
        sv, tv = sp.Rational(random.randint(2, 9), 11), sp.Rational(random.randint(2, 9), 13)
        rhs = 0
        for i in range(m):
            w = _prod([(X[i] - tv*X[j])/(X[i] - X[j]) for j in range(m) if j != i])
            rhs += X[i]*(1 + sv*X[i]*z)*w*Tk(k-1, X[:i] + X[i+1:], sv, tv)
        lhs = sum(tv**i for i in range(k))*Tk(k, X, sv, tv)
        good = sp.expand(lhs - rhs) == 0
        ok &= good
        print(f'm={m} k={k}', 'OK' if good else 'FAIL', flush=True)
print('ALL OK' if ok else 'SOME FAIL')
