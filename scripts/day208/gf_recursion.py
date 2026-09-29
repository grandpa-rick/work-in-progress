"""Day 208: symbolic outer-peel recursion (R)+(Lemma 2) for the two-variable GF
   Gamma_k(z,w) := sum_{a,b} z^a w^b t^{-C(k,2)} e_k(Y).(e_a e_b).
(R): [k] Gamma_k(X) = sum_i X_i(1+sX_i z)(1+sX_i w) prod_{j!=i} a_ij Gamma_{k-1}(X^_i),  Gamma_0 = E(z)E(w).
Representation: dict (lam, i, j) -> R(z,w,s,t)  meaning  R * e_lam * E(t^i z) E(t^j w).
Step: e_n(X^_i) = sum_c (-u)^c e_{n-c}, E(X^_i;g) = E(g)/(1+g u); Lemma 2: sum_i f(X_i) prod a_ij = Omega'[f]/(1-t)
for f(0)=0, Omega'[u^n] = q_n (q_0 = 1), and Omega'[1/(1+g u)] E(g) = Q(-g)E(g) = E(t g).
Only Q's that occur are single Q(-gamma) against its own E(gamma) -> absorbed as a chain shift.
"""
import sympy as sp, itertools, pickle, sys
s, t, z, w, u, y = sp.symbols('s t z w u y')
NMAX = 12
E = sp.symbols(f'e1:{NMAX+1}')
def e(n):
    return sp.Integer(1) if n == 0 else (E[n-1] if 0 < n <= NMAX else sp.Integer(0))
# q_n in e-basis
_Em = 1 + sum((-y)**n * E[n-1] for n in range(1, 8))
_Et = 1 + sum((-t*y)**n * E[n-1] for n in range(1, 8))
_Q = sp.expand(sp.series(_Et/_Em, y, 0, 8).removeO())
qn = [sp.expand(_Q.coeff(y, n)) for n in range(8)]

def lam_to_expr(lam):
    return sp.Mul(*[e(n) for n in lam])

def expr_to_terms(P):
    """polynomial in e's -> dict lam -> coeff"""
    P = sp.expand(P)
    out = {}
    for mon, c in sp.Poly(P, *E).terms():
        lam = []
        for idx, ex in enumerate(mon): lam += [idx+1]*ex
        lam = tuple(sorted(lam, reverse=True))
        out[lam] = out.get(lam, 0) + c
    return out

def step(G, k):
    new = {}
    def add(key, val):
        new[key] = new.get(key, 0) + val
    for (lam, i, j), R in G.items():
        g, d = t**i*z, t**j*w
        # e_lam(X^) as polynomial in u with e-coefficients
        ehat = sp.Integer(1)
        for n in lam:
            ehat *= sum((-u)**c * e(n-c) for c in range(n+1))
        ehat = sp.Poly(sp.expand(ehat), u)
        for (p,), cf in ehat.terms():
            N = sp.expand(u**(p+1)*(1+s*u*z)*(1+s*u*w))
            D = sp.expand((1+g*u)*(1+d*u))
            Pq, _ = sp.div(sp.Poly(N, u), sp.Poly(D, u))
            A = sp.cancel(N.subs(u, -1/g)/(1 - d/g))
            B = sp.cancel(N.subs(u, -1/d)/(1 - g/d))
            base = R*cf  # cf polynomial in e's
            # polynomial part: sum P_n q_n  (times E(g)E(d))
            polypart = sum(Pq.as_expr().coeff(u, n) * qn[n] for n in range(Pq.degree()+1)) if not Pq.is_zero else 0
            for mu, c2 in expr_to_terms(sp.expand(base*polypart)).items():
                add((mu, i, j), c2)
            for mu, c2 in expr_to_terms(sp.expand(base)).items():
                add((mu, i+1, j), c2*A)
                add((mu, i, j+1), c2*B)
    brk = sum(t**l for l in range(k))
    out = {}
    for key, val in new.items():
        v = sp.factor(sp.cancel(val/((1-t)*brk)))
        if v != 0: out[key] = v
    return out

if __name__ == '__main__':
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    G = {((), 0, 0): sp.Integer(1)}
    res = {0: G}
    for k in range(1, K+1):
        G = step(G, k)
        res[k] = G
        print(f'--- Gamma_{k}: {len(G)} terms', flush=True)
        for key in sorted(G, key=lambda x: (x[1], x[2], x[0])):
            print('  ', key, ':', G[key])
    pickle.dump(res, open(f'gf_Gamma_upto{K}.pkl', 'wb'))
