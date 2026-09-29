"""Folded human-readable form, a <= b:
  t^{-C(k,2)} e_k(Y).(e_a e_b) = sum_lam e_lam [ sum_{d1+d2=k-|lam|} Ch_{lam,d1,d2}(u,v) e_{a+d1} e_{b+d2}
                                              + sum_{q=0}^{min(b, a+k-|lam|-1)} T_lam(q) e_q e_{N_lam-q} ],
  u=t^a, v=t^b, U=t^q, N_lam=a+b+k-|lam|;  Ch from Laurent parts: Ch = sum_{i,j} u^i v^j t^{i d1 + j d2} L_{key,-d2}.
Prints formulas and verifies against flint AHA data."""
import sys, pickle, sympy as sp
from gf_recursion import s, t
U, u, v = sp.symbols('U u v')
k = int(sys.argv[1]); D = pickle.load(open(f'decomp_k{k}.pkl', 'rb')); T = pickle.load(open(f'tails_k{k}.pkl', 'rb'))
Ch = {}
for (lam, i, j), (L, P) in D.items():
    for n, c in L.items():
        d2 = -n; d1 = k - sum(lam) + n
        assert d2 >= 0 and d1 >= 0, (lam, i, j, n)
        Ch[(lam, d1, d2)] = Ch.get((lam, d1, d2), 0) + u**i * v**j * t**(i*d1 + j*d2) * c
Ch = {key: sp.factor(sp.cancel(c)) for key, c in Ch.items() if sp.cancel(c) != 0}
if __name__ == '__main__':
    print(f'=== k={k}: chain coefficients Ch[lam, d1, d2] (term e_lam e_(a+d1) e_(b+d2)) ===')
    for key in sorted(Ch): print('  ', key, ':', Ch[key])
    print(f'=== k={k}: tail densities T_lam(q) (term e_lam e_q e_(N-q)), U=t^q ===')
    for lam in sorted(T): print('  ', lam, ':', T[lam])
    if len(sys.argv) > 2:
        data = pickle.load(open(sys.argv[2], 'rb')); ok = True
        for (a, b), d in sorted(data.items()):
            sub = {u: t**a, v: t**b}
            out = {}
            def add(parts, c):
                mu = tuple(sorted([m for m in parts if m > 0], reverse=True))
                if min(parts) < 0: return
                out[mu] = out.get(mu, 0) + c
            for (lam, d1, d2), c in Ch.items(): add(lam + (a+d1, b+d2), c.subs(sub))
            for lam, Tl in T.items():
                N = a + b + k - sum(lam)
                for q in range(0, min(b, a + k - sum(lam) - 1) + 1):
                    add(lam + (q, N-q), Tl.subs(sub).subs(U, t**q))
            good = all(sp.cancel(out.get(l, 0) - d.get(l, 0)) == 0 for l in set(out) | set(d))
            ok &= good
            print(f'  folded formula vs AHA (a,b)=({a},{b}): {"OK" if good else "FAIL"}', flush=True)
        print('ALL OK' if ok else 'SOME FAIL')
