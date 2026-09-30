"""Day 210: ell-column generalisation of scripts/day208/gf_recursion.py (outer-peel recursion (R_ell) + Lemma 2).
  Gamma_k(z_1..z_ell) := sum_a z^a t^{-C(k,2)} e_k(Y).(e_{a_1}...e_{a_ell}).
(R_ell): [k] Gamma_k(X) = sum_i X_i prod_c(1+sX_i z_c) prod_{j!=i} a_ij Gamma_{k-1}(X^_i),  Gamma_0 = prod_c E(z_c).
Representation: dict (lam, (i_1..i_ell)) -> R  meaning  R * e_lam * prod_c E(t^{i_c} z_c).
Same step as day208 (verbatim logic), with ell denominators (1+gam_c u) instead of two.
"""
import sympy as sp, pickle, sys
s, t, u, y = sp.symbols('s t u y')
NMAX = 12
E = sp.symbols(f'e1:{NMAX+1}')
def e(n):
    return sp.Integer(1) if n == 0 else (E[n-1] if 0 < n <= NMAX else sp.Integer(0))
_Em = 1 + sum((-y)**n * E[n-1] for n in range(1, 9))
_Et = 1 + sum((-t*y)**n * E[n-1] for n in range(1, 9))
_Q = sp.expand(sp.series(_Et/_Em, y, 0, 9).removeO())
qn = [sp.expand(_Q.coeff(y, n)) for n in range(9)]

def expr_to_terms(P):
    P = sp.expand(P)
    out = {}
    if P == 0: return out
    for mon, c in sp.Poly(P, *E).terms():
        lam = []
        for idx, ex in enumerate(mon): lam += [idx+1]*ex
        lam = tuple(sorted(lam, reverse=True))
        out[lam] = out.get(lam, 0) + c
    return out

def step(G, k, Z):
    ell = len(Z)
    new = {}
    def add(key, val): new[key] = new.get(key, 0) + val
    for (lam, I), R in G.items():
        gam = [t**I[c]*Z[c] for c in range(ell)]
        ehat = sp.Integer(1)
        for n in lam:
            ehat *= sum((-u)**c * e(n-c) for c in range(n+1))
        ehat = sp.Poly(sp.expand(ehat), u)
        D = sp.expand(sp.Mul(*[1+g*u for g in gam]))
        for (p,), cf in ehat.terms():
            N = sp.expand(u**(p+1)*sp.Mul(*[1+s*u*zc for zc in Z]))
            Pq, _ = sp.div(sp.Poly(N, u), sp.Poly(D, u))
            res = [sp.cancel(N.subs(u, -1/gam[c]) / sp.Mul(*[1-gam[d]/gam[c] for d in range(ell) if d != c])) for c in range(ell)]
            base = R*cf
            polypart = sum(Pq.as_expr().coeff(u, n)*qn[n] for n in range(Pq.degree()+1)) if not Pq.is_zero else 0
            for mu, c2 in expr_to_terms(base*polypart).items():
                add((mu, I), c2)
            for mu, c2 in expr_to_terms(base).items():
                for c in range(ell):
                    J = list(I); J[c] += 1
                    add((mu, tuple(J)), c2*res[c])
    brk = sum(t**l for l in range(k))
    out = {}
    for key, val in new.items():
        v = sp.factor(sp.cancel(val/((1-t)*brk)))
        if v != 0: out[key] = v
    return out

def run(ell, K):
    Z = sp.symbols(f'z1:{ell+1}')
    G = {((), (0,)*ell): sp.Integer(1)}
    res = {0: G}
    for k in range(1, K+1):
        G = step(G, k, Z); res[k] = G
        print(f'--- ell={ell} Gamma_{k}: {len(G)} terms; lam lengths present: {sorted({len(l) for l, _ in G})}', flush=True)
    return Z, res

if __name__ == '__main__':
    ell = int(sys.argv[1]); K = int(sys.argv[2])
    Z, res = run(ell, K)
    pickle.dump(res, open(f'/home/agent/projects/scripts/day210/gamma_ell{ell}_k{K}.pkl', 'wb'))
