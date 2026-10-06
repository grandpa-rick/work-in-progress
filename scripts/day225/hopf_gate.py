# Day 225 Hopf-route gate. grade: computed.
# <H,p_x> = H(1,z,...,z^{x-1}), z = prim. x-th root of unity, for H homogeneous deg x (only p_(x) survives).
# Phi_a(f;x,y) = coefficient of u^x in (T_a f)(u*Zx ⊔ Zy), T_a g = sum_{|A|=a} c_A X_A g(X_A).
# Piece i = sum over A with |A∩X|=i. Pieces are rational in u (poles at |u|=1): we take
#   'in'  = Taylor coeff at u=0 (|X|<<|Y|), 'out' = Laurent coeff at u=inf (|X|>>|Y|), via FFT on circles r=0.5, 2.
import numpy as np, itertools, sys, sympy
from fractions import Fraction as Fr
t_s = sympy.symbols('t')
def fvals(name, xs):
    xs = np.array(xs)
    if name == '1': return 1.0
    if name == 'e1': return xs.sum()
    if name == 'p2': return (xs**2).sum()
    if name == 'e1^2': return xs.sum()**2
    if name == 'e2': return sum(xs[i]*xs[j] for i in range(len(xs)) for j in range(i+1, len(xs)))
    if name == 'e3': return sum(xs[i]*xs[j]*xs[k] for i,j,k in itertools.combinations(range(len(xs)),3))
    if name == 'p3': return (xs**3).sum()
    if name == 'e1e2': return fvals('e1',xs)*fvals('e2',xs)
DEG = {'1':0,'e1':1,'p2':2,'e1^2':2,'e2':2,'e3':3,'p3':3,'e1e2':3}
def pieces_at(a, f, x, y, t, u):
    X = [u*np.exp(2j*np.pi*k/x) for k in range(x)]; Y = [np.exp(2j*np.pi*k/y) for k in range(y)]
    V = X + Y; N = len(V); out = np.zeros(a+1, complex)
    for A in itertools.combinations(range(N), a):
        S = set(A); c = 1.0
        for i in A:
            for j in range(N):
                if j not in S: c *= (V[i]-t*V[j])/(V[i]-V[j])
        xa = [V[i] for i in A]
        out[sum(1 for i in A if i < x)] += c*np.prod(xa)*fvals(f, xa)
    return out
def coeffs(a, f, x, y, t, M=48):
    res = {}
    for tag, r in (('in', 0.5), ('out', 2.0)):
        us = r*np.exp(2j*np.pi*np.arange(M)/M)
        vals = np.array([pieces_at(a, f, x, y, t, u) for u in us])  # M x (a+1)
        c = (vals * (us[:, None]**(-x))).mean(axis=0)
        res[tag] = c
    return res
def fit(fun, D, _):
    # t on K-th roots of unity -> exact poly coefficients by DFT; ok = top coeffs vanish (no aliasing)
    K = 40
    ts_ = np.exp(2j*np.pi*np.arange(K)/K)
    ys = np.array([fun(tt) for tt in ts_])
    c = np.fft.fft(ys)/K   # c[m] = sum_k ys_k w^{-km}/K  -> coefficient of t^m
    ok = np.all(np.abs(c[D+1:]) < 1e-6)
    cs = [Fr(float(v.real)).limit_denominator(2000) for v in c[:D+1]]
    ok = ok and all(abs(v.imag) < 1e-6 for v in c) and all(abs(float(q) - v.real) < 1e-6 for q, v in zip(cs, c[:D+1]))
    poly = sum(sympy.Rational(q.numerator, q.denominator)*t_s**m for m, q in enumerate(cs))
    return sympy.factor(poly), ok
if __name__ == '__main__':
    cases = eval(sys.argv[1])
    for a, f, x, y in cases:
        n = a + DEG[f]; assert x + y == n
        D = 2*n + 8
        cache = {}
        def get(tt):
            if tt not in cache: cache[tt] = coeffs(a, f, x, y, tt)
            return cache[tt]
        norm = lambda tt: (1-tt**x)*(1-tt**y)
        get = lambda tt, g=get: g(complex(tt))
        tot, ok = fit(lambda tt: get(tt)['in'].sum()*norm(tt), D, None)
        print(f'a={a} f={f} (x,y)=({x},{y})  TOTAL*(1-t^x)(1-t^y) = {tot}  [fit ok={ok}]', flush=True)
        for tag in ('in', 'out'):
            for i in range(a+1):
                p, ok = fit(lambda tt: get(tt)[tag][i]*norm(tt), D, None)
                print(f'    {tag:3s} i=|A_X|={i}: piece*(1-t^x)(1-t^y) = {p}  [ok={ok}]', flush=True)
