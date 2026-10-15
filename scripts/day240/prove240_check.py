"""PROVE 240: reconcile Y^H_i = theta^{N-1}(Y_{i,-1})^{-1} (same q,theta) with rem:N-DFK (E_k = t^{k(m-k)} M_{k;1} at (q,t)_DFK=(s,t^-1)).
Printed sources: DFK 1704.00154 (1.2) Gamma_i: x_i -> q x_i; (1.5) M_{a;n}=sum_I x_I^n prod_{i in I, j notin I}(t x_i-x_j)/(x_i-x_j) Gamma_I;
Def D_{a;n}=q^{-an} e_a(Y_{.,n}); Thm genmacthm: D_{a;n}|_S = theta^{-a(N-a)} M_{a;n}.
Long version: a_ij=(x_i-t x_j)/(x_i-x_j), c_A, E_k F = sum_A c_A x_A F|_{x_i->s x_i, i in A}; thm:subset e_k(Y^H)=t^{C(k,2)} E_k; s=q_H^{-1} (l.523).
"""
import sympy as sp, itertools
from clio_ymatch import th, q, mk, simp, H_Y, D_Yn, D_Yn_inv, prod
t = th**2
s = sp.Symbol('s'); T = sp.Symbol('T'); Q = sp.Symbol('Q')   # T,Q: DFK's t,q as free symbols for printed formulas

def symbasis(x, deg):
    N = len(x); out = []
    for lam in itertools.product(range(deg+1), repeat=N):
        if list(lam) == sorted(lam, reverse=True) and 0 < sum(lam) <= deg:
            out.append(sp.expand(sum(prod([v**e for v, e in zip(p, lam)]) for p in set(itertools.permutations(x)))))
    return out

def E_k(F, x, k, sv):
    N = len(x); r = 0
    for A in itertools.combinations(range(N), k):
        c = prod([(x[i]-t*x[j])/(x[i]-x[j]) for i in A for j in range(N) if j not in A])
        r += c*prod([x[i] for i in A])*F.subs({x[i]: sv*x[i] for i in A}, simultaneous=True)
    return sp.cancel(r)

def M_printed(F, x, a, n, qq, tt):   # DFK (1.5) with (1.2)
    N = len(x); r = 0
    for I in itertools.combinations(range(N), a):
        c = prod([(tt*x[i]-x[j])/(x[i]-x[j]) for i in I for j in range(N) if j not in I])
        r += prod([x[i] for i in I])**n*c*F.subs({x[i]: qq*x[i] for i in I}, simultaneous=True)
    return sp.cancel(r)

def ek_ops(F, x, k, op):   # e_k of commuting operators op(f,i) applied to F
    N = len(x); r = 0
    for I in itertools.combinations(range(1, N+1), k):
        g = F
        for i in I: g = op(g, i)
        r += g
    return sp.cancel(r)

def z(e): return simp(e) == 0

if __name__ == "__main__":
  pass
for N, deg in ([] if __name__ != "__main__" else [(2, 3), (3, 3), (4, 2)]):
    x = mk(N); Fs = symbasis(x, deg)
    print(f"N={N}: {len(Fs)} monomial symmetric fns, deg<={deg}")
    for k in range(1, N+1):
        # T2: thm:subset with Hikita's operators, s=q^-1
        b2 = sum(not z(ek_ops(F, x, k, lambda g, i: H_Y(g, x, i, t)) - t**sp.binomial(k, 2)*E_k(F, x, k, 1/q)) for F in Fs)
        # T3: printed rem:N-DFK: E_k = t^{k(m-k)} M_{k;1} at (q,t)_DFK=(s,t^-1)
        b3 = sum(not z(E_k(F, x, k, s) - t**(k*(N-k))*M_printed(F, x, k, 1, s, 1/t)) for F in Fs)
        # T3': DFK operators themselves at (q,theta)_DFK=(s, theta^-1): q^{-k} e_k(Y_{.,1}) = theta_D^{-k(N-k)} M_{k;1}
        def Y1inv(g, i):  # DFK Y_{i,1} at q->s, theta->1/theta
            return D_Yn(g, x, i, 1).subs({q: s, th: 1/th}, simultaneous=True)
        b3p = sum(not z(s**(-k)*ek_ops(F, x, k, Y1inv) - th**(k*(N-k))*M_printed(F, x, k, 1, s, 1/t)) for F in Fs)
        # T4: bridge (I): M_{N;-1}^{-1} M_{N-k;-1} at (Q,T) = Q^{-k} T^{k(N-k)} M_{k;1} at (Q^-1,T^-1)
        def MNinv(F): return sp.cancel(Q**(-N)*prod(x)*F.subs({v: v/Q for v in x}, simultaneous=True))
        b4 = sum(not z(MNinv(M_printed(F, x, N-k, -1, Q, T)) - Q**(-k)*T**(k*(N-k))*M_printed(F, x, k, 1, 1/Q, 1/T)) for F in Fs)
        # controls for (I): invert only q / only t / no k<->N-k
        c1 = sum(not z(MNinv(M_printed(F, x, N-k, -1, Q, T)) - Q**(-k)*T**(k*(N-k))*M_printed(F, x, k, 1, 1/Q, T)) for F in Fs)
        c2 = sum(not z(MNinv(M_printed(F, x, N-k, -1, Q, T)) - Q**(-k)*T**(k*(N-k))*M_printed(F, x, k, 1, Q, 1/T)) for F in Fs)
        c3 = sum(not z(MNinv(M_printed(F, x, k, -1, Q, T)) - Q**(-k)*T**(k*(N-k))*M_printed(F, x, k, 1, 1/Q, 1/T)) for F in Fs) if 2*k != N else 'n/a'
        # T5: rem:N-DFK route vs Hikita directly: e_k(Y^H) = t^{C(k,2)} t^{k(N-k)} M_{k;1}(s=q^-1, t^-1)
        b5 = sum(not z(ek_ops(F, x, k, lambda g, i: H_Y(g, x, i, t)) - t**(sp.binomial(k, 2)+k*(N-k))*M_printed(F, x, k, 1, 1/q, 1/t)) for F in Fs)
        print(f"  k={k}: T2 subset bad={b2}  T3 printed E_k=t^(k(m-k))M_k;1 bad={b3}  T3' DFK-ops bad={b3p}  T4 bridge(I) bad={b4}  "
              f"controls(I) q-only={c1} t-only={c2} no-k-swap={c3}  T5 bad={b5}")
    # T6: operator-level n=+1 at inverted parameters? Y^H_i o Y^DFK_{j,1}|_{(q^-1,theta^-1)} and Y^H_i o (Y_{j,1})^-1|_{inv}: multiplication op?
    ms = [prod([v**e for v, e in zip(x, ex)]) for ex in itertools.product(range(3), repeat=N) if sum(ex) <= 2]
    inv = {q: 1/q, th: 1/th}
    for i in range(1, N+1):
        for j in range(1, N+1):
            for name, op in [('Y', lambda f: D_Yn(f, x, j, 1)), ('Yinv', lambda f: D_Yn_inv(f, x, j, 1))]:
                comp = lambda f: H_Y(sp.cancel(op(f).subs(inv, simultaneous=True)), x, i, t)
                g = simp(comp(sp.Integer(1)))
                bad = sum(not z(comp(f) - g*f) for f in ms)
                if bad == 0: print(f"  T6 HIT: Y^H_{i} o {name}_({j},1)|inv is mult by {g}")
    print(f"  T6 done (only HIT lines would indicate an operator-level n=+1 match at inverted params)")
