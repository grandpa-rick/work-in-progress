"""Fast closed-form coefficient extraction. For each key, apart(rho(x)) = sum_{d} L_d x^d + sum_c A_c/(x - x_c), so
[x^n] rho = L_n + [n>=0] * sum_c (-A_c) x_c^{-n-1}.  x_c in {1, t, 1/t, ...}.
Formula:  t^{-C(k,2)} e_k(Y).(e_a e_b) = sum_{(lam,i,j)} sum_{p+q=a+b+k-|lam|} t^{ip+jq} rho_key[b-q] e_lam e_p e_q."""
import sys, pickle, sympy as sp
from gf_recursion import s, t, z, w
x = sp.Symbol('x'); nn = sp.Symbol('n')
def decompose(R):
    r = sp.cancel(R.subs({z: 1, w: x}))
    L = {}; P = []
    for term in sp.Add.make_args(sp.apart(r, x)):
        num, den = sp.fraction(sp.factor(term))
        dp = sp.Poly(den, x)
        if dp.degree() == 0 or (len(dp.terms()) == 1):  # monomial denominator -> Laurent term
            tt = sp.expand(term)
            for mono in sp.Add.make_args(tt):
                c, d = mono.as_coeff_exponent(x) if mono.has(x) else (mono, 0)
                if mono.has(x):
                    d = sp.degree(sp.fraction(mono)[0], x) - sp.degree(sp.fraction(mono)[1], x)
                    c = sp.cancel(mono / x**d)
                L[int(d)] = L.get(int(d), 0) + c
        else:
            assert dp.degree() == 1, term
            x0 = sp.solve(den, x)[0]
            A = sp.cancel(term * (x - x0))
            P.append((sp.simplify(A), sp.simplify(x0)))
    return L, P
def coeff_fun(L, P, n):
    v = L.get(n, 0)
    if n >= 0: v += sum(-A * x0**(-n-1) for A, x0 in P)
    return v
def extract(D, k, a, b):
    out = {}
    for (lam, i, j), (L, P) in D.items():
        tot = a + b + k - sum(lam)
        for q in range(0, tot+1):
            p = tot - q
            c = coeff_fun(L, P, b - q)
            if c == 0: continue
            mu = tuple(sorted([m for m in lam + (p, q) if m > 0], reverse=True))
            out[mu] = out.get(mu, 0) + t**(i*p + j*q) * c
    return out
if __name__ == '__main__':
    k = int(sys.argv[1]); G = pickle.load(open(sys.argv[2], 'rb'))[k]; data = pickle.load(open(sys.argv[3], 'rb'))
    D = {key: decompose(R) for key, R in G.items()}
    pickle.dump(D, open(f'decomp_k{k}.pkl', 'wb'))
    ok = True
    for (a, b), d in sorted(data.items()):
        ex = extract(D, k, a, b)
        good = all(sp.cancel(ex.get(l, 0) - d.get(l, 0)) == 0 for l in set(ex) | set(d))
        ok &= good
        print(f'k={k} (a,b)=({a},{b}): {"OK" if good else "FAIL"}', flush=True)
    print('ALL OK' if ok else 'SOME FAIL')
