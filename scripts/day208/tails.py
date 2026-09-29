"""Tail densities: T_lam(q) = sum_{keys (lam,i,j)} t^{ip+jq} G_key(b-q)  (pole part, valid for q<=b),
in U=t^q, u=t^a, v=t^b; N_lam = a+b+k-|lam|, t^p = t^{N_lam}/U. Also antisymmetry test T(q)+T(N-q)=0."""
import sys, pickle, sympy as sp
from gf_recursion import s, t
U, u, v, V = sp.symbols('U u v V')
k = int(sys.argv[1]); D = pickle.load(open(f'decomp_k{k}.pkl', 'rb'))
lams = sorted({key[0] for key in D})
out = {}
for lam in lams:
    Tq = 0
    tN = u*v*t**(k - sum(lam))
    for (l2, i, j), (L, P) in D.items():
        if l2 != lam: continue
        # G(n) = sum -A x0^{-n-1}, n=b-q: x0^{-n-1} = x0^{-1} * x0^{-b} * x0^{q}
        for A, x0 in P:
            e = [ee for ee in range(-4, 5) if sp.simplify(x0 - t**ee) == 0][0]  # x0 = t^e
            Tq += -A * t**(-e) * v**(-e) * U**e * (tN/U)**i * U**j
    Tq = sp.factor(sp.cancel(Tq))
    anti = sp.cancel(Tq + Tq.subs(U, tN/U, simultaneous=True))
    print('lam =', lam, ' T(q) =', Tq, '   T(q)+T(N-q) =', anti)
    out[lam] = Tq
pickle.dump(out, open(f'tails_k{k}.pkl', 'wb'))
