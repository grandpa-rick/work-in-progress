# Day 236 print check: re-implement Thm 6.2 (N), Prop 6.7 (NS), Thm A.4/Cor A.5 from the PRINTED text
# of longversion.pdf (sec 6 + App A). Independent of scripts/day216c engines.
import itertools, sys
from sympy import symbols, Rational as R, expand, cancel, Poly, Matrix, binomial, factor, simplify, S
x = symbols('x1:6'); s_, t_ = symbols('s t')

def sub_swap(f, i):  # s_i, 1-based i swaps x_i,x_{i+1}
    a, b = x[i-1], x[i]
    return f.xreplace({a: b, b: a})

def T(f, i, t):
    sf = sub_swap(f, i)
    return expand(t*sf + (t-1)*cancel(x[i]*(sf - f)/(x[i-1]-x[i])))

def Tinv(f, i, t):
    return expand(T(f, i, t)/t - (1-1/t)*f)

def pi0(f, m, s):  # f(x2..xm, s x1)
    rep = {x[j]: x[j+1] for j in range(m-1)}; rep[x[m-1]] = s*x[0]
    return expand(f.xreplace(rep))

def pitw(f, m, s): return expand(x[0]*pi0(f, m, s))

def Yop(f, i, m, s, t, pi):  # Y_i = t^{m-i} T_{i-1}..T_1 pi T_{m-1}^{-1}..T_i^{-1}; rightmost first
    g = f
    for j in range(i, m): g = Tinv(g, j, t)        # T_i^{-1} first, ... T_{m-1}^{-1}
    g = pi(g, m, s)
    for j in range(1, i): g = T(g, j, t)  # T_1 acts first, T_{i-1} last
    return expand(t**(m-i)*g)

def comps(d, m): return [c for c in itertools.product(range(d+1), repeat=m) if sum(c) == d]
def mono(c, m):
    r = S(1)
    for j in range(m): r *= x[j]**c[j]
    return r
def coeffs(f, m, d):
    P = Poly(f, *x[:m]); return {c: P.coeff_monomial(mono(c, m)) for c in comps(d, m)}

# ---------- Part C: printed example and triangularity (symbolic s,t) ----------
def partC():
    m = 3; ok = True
    lam = (2, 0, 1)
    got = Yop(mono(lam, m), 1, m, s_, t_, pi0)
    printed = s_**2*mono((2,0,1),m) - s_*(t_-1)*mono((1,0,2),m) - t_*(t_-1)*mono((0,2,1),m) - s_*(t_-1)*mono((1,1,1),m)
    e = expand(got - printed) == 0; print('C1 printed example Y1o x^(2,0,1):', e); ok &= e
    def sortp(c): return tuple(sorted(c, reverse=True))
    def dom_le(a, b):
        A = B = 0
        for u, v in zip(a, b):
            A += u; B += v
            if A > B: return False
        return True
    def prec(mu, la):
        if sortp(mu) != sortp(la): return dom_le(sortp(mu), sortp(la))
        return mu != la and dom_le(mu, la)
    def rk(la, i): return sum(1 for j in range(m) if la[j] > la[i]) + sum(1 for j in range(i) if la[j] == la[i])
    bad = 0; n = 0
    for d in range(4):
        for la in comps(d, m):
            for i in range(1, m+1):
                C = coeffs(Yop(mono(la, m), i, m, s_, t_, pi0), m, d)
                for mu, c in C.items():
                    c = expand(c)
                    if mu == la:
                        if expand(c - s_**la[i-1]*t_**rk(la, i-1)) != 0: bad += 1
                    elif c != 0:
                        n += 1
                        if not prec(mu, la): bad += 1
    print(f'C2 Thm A.4 m=3 d<=3: diag=s^lam_i t^rk_i and off-diag prec lam: bad={bad} (off-diag terms {n})'); ok &= bad == 0
    return ok

# ---------- Part A: Thm 6.2 at rational (s,t) ----------
def partitions(n, maxlen):
    def gen(n, mx):
        if n == 0: yield (); return
        for k in range(min(n, mx), 0, -1):
            for r in gen(n-k, k): yield (k,)+r
    return [p for p in gen(n, n) if len(p) <= maxlen]
def msym(nu, m):
    nu = tuple(nu) + (0,)*(m-len(nu))
    return sum(mono(c, m) for c in set(itertools.permutations(nu)))
def nfun(nu): return sum(i*v for i, v in enumerate(nu))
def conj(nu): return tuple(sum(1 for v in nu if v > j) for j in range(nu[0])) if nu else ()
def dom_le_p(a, b):
    A = B = 0
    for j in range(max(len(a), len(b))):
        A += a[j] if j < len(a) else 0; B += b[j] if j < len(b) else 0
        if A > B: return False
    return True

def macdonaldP(n, m, q, tau):
    # Macdonald VI (3.4): D = sum_i prod_{j!=i} (tau x_i - x_j)/(x_i - x_j) T_{q,x_i}; P_nu = m_nu + lower, eigenvector.
    parts = partitions(n, m); basis = [msym(p, m) for p in parts]
    def D(f):
        r = 0
        for i in range(m):
            A = S(1)
            for j in range(m):
                if j != i: A *= (tau*x[i]-x[j])/(x[i]-x[j])
            r += A*f.xreplace({x[i]: q*x[i]})
        return expand(cancel(together_(r)))
    def lead_coeffs(f):  # coefficients in m-basis via leading monomials
        P = Poly(f, *x[:m]); return [P.coeff_monomial(mono(tuple(p)+(0,)*(m-len(p)), m)) for p in parts]
    M = Matrix([lead_coeffs(D(b)) for b in basis]).T  # column j = D(m_j) in m-basis
    out = {}
    for j, nu in enumerate(parts):
        ev = sum(q**(nu[i] if i < len(nu) else 0)*tau**(m-1-i) for i in range(m))
        # solve (M - ev) v = 0 with v_nu = 1, support on mu <= nu
        idx = [k for k, mu in enumerate(parts) if dom_le_p(mu, nu)]
        A = (M - ev*Matrix.eye(len(parts)))
        others = [k for k in idx if k != j]
        if others:
            sub = A.extract(list(range(len(parts))), others); rhs = -A.extract(list(range(len(parts))), [j])
            sol = sub.solve_least_squares(rhs) if sub.rows != sub.cols else sub.solve(rhs)
            assert (sub*sol - rhs).is_zero_matrix
            v = {k: sol[a] for a, k in enumerate(others)}
        else: v = {}
        v[j] = 1
        assert (A*Matrix([v.get(k, 0) for k in range(len(parts))])).is_zero_matrix
        out[nu] = expand(sum(v[k]*basis[k] for k in v))
    return out, parts, lead_coeffs
from sympy import together as together_

def partA(m, sv, tv, maxtot, neg=False):
    ok = True; q, tau = sv, 1/tv
    Ps = {}; LC = {}
    for n in range(maxtot+1):
        Ps[n], parts, lc = macdonaldP(n, m, q, tau); LC[n] = (parts, lc)
    def toP(f, n):  # expand symmetric f of degree n in P basis (triangular)
        parts, lc = LC[n]; Mx = Matrix([lc(Ps[n][p]) for p in parts]).T
        c = Mx.solve(Matrix(lc(f))); return {p: c[a] for a, p in enumerate(parts)}
    def Nop(f, n, inv=False):
        c = toP(f, n); w = (lambda p: tv**nfun(p)*sv**nfun(conj(p))) if not neg else (lambda p: tv**nfun(conj(p))*sv**nfun(p))
        return expand(sum(c[p]*(w(p)**(-1 if inv else 1))*Ps[n][p] for p in c))
    def Ek(f, k):
        r = 0
        for A in itertools.combinations(range(m), k):
            cA = S(1)
            for i in A:
                for j in range(m):
                    if j not in A: cA *= (x[i]-tv*x[j])/(x[i]-x[j])
            xA = S(1)
            for i in A: xA *= x[i]
            r += cA*xA*f.xreplace({x[i]: sv*x[i] for i in A})
        return expand(cancel(together_(r)))
    def ek(k): return sum(mono(tuple(1 if j in A else 0 for j in range(m)), m) for A in itertools.combinations(range(m), k))
    cnt = 0; bad = 0
    for n in range(maxtot):
        for nu in Ps[n]:
            for k in range(1, min(m, maxtot-n)+1):
                F = Ps[n][nu]
                lhs = Ek(F, k)
                rhs = expand(tv**(-binomial(k, 2))*Nop(expand(ek(k)*Nop(F, n, inv=True)), n+k))
                cnt += 1
                if expand(lhs-rhs) != 0:
                    bad += 1
                    if not neg: print('  FAIL', nu, k)
    # negative control: q = 1/s instead of s in P (should fail somewhere)
    print(f'A Thm 6.2{" NEG(n(nu),n(nu^prime) swapped)" if neg else ""} m={m} (s,t)=({sv},{tv}) |nu|+k<={maxtot}: {cnt-bad}/{cnt} OK'); ok &= bad == 0
    return ok

# ---------- Part B: Prop 6.7 at rational point, m=3, d<=2 ----------
def partB(m, sv, tv, dmax, thetaspec=(1, 1)):
    sa, sb = thetaspec  # theta(s^a t^r) = s^{sa*C(a,2)} t^{sb*a r}
    def Ymat(i, d, pi):
        cs = comps(d, m); return cs, Matrix([[coeffs(Yop(mono(c, m), i, m, sv, tv, pi), m, d)[c2] for c in cs] for c2 in cs])
    # eigenbasis of Y°: Z = sum c^{i-1} Y°_i
    def eig(d):
        cs, _ = Ymat(1, d, pi0)
        Ys = [Ymat(i, d, pi0)[1] for i in range(1, m+1)]
        Z = Ys[0]
        for i in range(1, m): Z = Z + R(3)**i*Ys[i]
        def rk(la, i): return sum(1 for j in range(m) if la[j] > la[i]) + sum(1 for j in range(i) if la[j] == la[i])
        vecs = []; gam = []
        for la in cs:
            y = [R(sv)**la[i]*R(tv)**rk(la, i) for i in range(m)]
            z = sum(R(3)**i*y[i] for i in range(m))
            ns = (Z - z*Matrix.eye(len(cs))).nullspace(); assert len(ns) == 1, (la, len(ns))
            v = ns[0]
            for i in range(m): assert (Ys[i]*v - y[i]*v).is_zero_matrix
            vecs.append(v)
            g = 1
            for i in range(m):
                a, r = la[i], rk(la, i); g *= R(sv)**(sa*binomial(a, 2))*R(tv)**(sb*a*r)
            gam.append(g)
        Emat = Matrix.hstack(*vecs); G = Emat*Matrix.diag(*gam)*Emat.inv()
        return cs, G
    ok = True; Gs = {d: eig(d) for d in range(dmax+2)}
    for d in range(dmax+1):
        csd, Gd = Gs[d]; cse, Ge = Gs[d+1]
        for i in range(1, m+1):
            Xi = Matrix([[coeffs(expand(x[i-1]*mono(c, m)), m, d+1)[c2] for c in csd] for c2 in cse])
            Yi = Matrix([[coeffs(Yop(mono(c, m), i, m, sv, tv, pitw), m, d+1)[c2] for c in csd] for c2 in cse])
            good = (Ge*Xi*Gd.inv() - Yi).is_zero_matrix
            ok &= good
    return ok

if __name__ == '__main__':
    allok = True
    r = partC(); allok &= r
    r = partB(3, R(5,7), R(-3,11), 2); print('B Prop 6.7 m=3 d<=2 (5/7,-3/11) gamma x_i gamma^-1 = Y_i (twisted), all i:', r); allok &= r
    r = partB(2, R(2,5), R(7,3), 3); print('B Prop 6.7 m=2 d<=3 (2/5,7/3):', r); allok &= r
    for spec in [(-1, 1), (1, -1), (-1, -1)]:
        r = partB(2, R(2,5), R(7,3), 1, spec); print(f'B NEG wrong Gaussian {spec} m=2 d<=1 (must be False):', r); allok &= (not r)
    r = partA(3, R(5,11), R(7,3), 5); allok &= r
    r = partA(4, R(3,7), R(-5,2), 4); allok &= r
    r = partA(3, R(5,11), R(7,3), 3, neg=True); allok &= (not r)
    print('ALL:', allok)
