"""Coefficient extraction from Gamma_k:  R_key(z,w) is homogeneous of degree |lam|-k, so R = z^{|lam|-k} rho(w/z),
rho(x) = R(1,x).  Then
   t^{-C(k,2)} e_k(Y).(e_a e_b) = sum_{key=(lam,i,j)} sum_{p+q=a+b+k-|lam|} t^{ip+jq} [x^{b-q}]rho_key(x) * e_lam e_p e_q
([x^n] = Laurent coefficient at x=0, i.e. region |w|<|z|).  Compare with flint AHA data."""
import sys, pickle, sympy as sp
from gf_recursion import s, t, z, w
x = sp.Symbol('x')
def rho(R): return sp.cancel(R.subs({z: 1, w: x}))
def laurent_coeff(r, n, cache={}):
    key = (r, n)
    if key not in cache:
        num, den = sp.fraction(sp.cancel(r))
        # den = x^o * d0(x), d0(0) != 0
        o = 0
        dp = sp.Poly(den, x)
        while dp.eval(0) == 0:
            dp = sp.Poly(sp.cancel(dp.as_expr()/x), x); o += 1
        N = n + o
        if N < 0: cache[key] = sp.Integer(0)
        else:
            ser = sp.series(num/dp.as_expr(), x, 0, N+1).removeO()
            cache[key] = sp.expand(ser).coeff(x, N)
    return cache[key]
def extract(G, k, a, b):
    out = {}
    for (lam, i, j), R in G.items():
        r = rho(R)
        tot = a + b + k - sum(lam)
        for q in range(0, tot+1):
            p = tot - q
            c = laurent_coeff(r, b - q)
            if c == 0: continue
            mu = tuple(sorted([n for n in lam + (p, q) if n > 0], reverse=True))
            out[mu] = out.get(mu, 0) + t**(i*p + j*q) * c
    return {m_: sp.factor(sp.simplify(v)) for m_, v in out.items() if sp.simplify(v) != 0}
if __name__ == '__main__':
    k = int(sys.argv[1]); G = pickle.load(open(sys.argv[2], 'rb'))[k]; data = pickle.load(open(sys.argv[3], 'rb'))
    s2 = sp.Symbol('s'); ok = True
    for (a, b), d in sorted(data.items()):
        ex = extract(G, k, a, b)
        good = all(sp.simplify(ex.get(l, 0) - d.get(l, 0)) == 0 for l in set(ex) | set(d))
        ok &= good
        print(f'k={k} (a,b)=({a},{b}): {"OK" if good else "FAIL"}', flush=True)
    print('ALL OK' if ok else 'SOME FAIL')
