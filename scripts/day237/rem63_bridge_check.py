# Day 237 (a): print check of longversion Rem 6.3 (rem:N-nabla), taken FROM THE PDF TEXT.
# Independent of scripts/day216*: symmetric functions as polynomials in power sums p1..pn,
#   P_nu  : Gram-Schmidt in the q,t-Hall pairing, dominance (total for n<=5) [Mac VI (4.7)]
#   J_nu  : c_nu(q,t) P_nu, c_nu = prod_s (1 - q^a t^{l+1})              [Mac VI (8.3)]
#   Ht_nu : Garsia-Haiman characterisation (NOT the printed normalisation), then the printed
#           normalisation Ht = t^{n(nu)} J[X/(1-1/t); q, 1/t] is CHECKED against it.
# Printed claims tested (pdftotext longversion.pdf, eb12457):
#   (P1) Ht_nu(x;q,t) = t^{n(nu)} J_nu[X/(1-t^{-1}); q, t^{-1}]
#   (P2) phi^{-1} Ht_nu(x;s,t) = (-1)^{|nu|} t^{|nu|+n(nu)} J_nu(x;s,t^{-1}),  phi(F)=F[-eps X/(1-t)]
#   (P3) N = phi^{-1} nabla_{q=s} phi  on Lambda^n,  N P_nu = t^{n(nu)} s^{n(nu')} P_nu, P_nu = P_nu(x;s,1/t)
#   (P4) t^{n(lam')} e*_lam = phi^{-1} nabla_{q=s} phi(e_lam), with e*_lam built from the printed
#        Sec.2 operators T_i, pi, Y_i and (eq:star-reading) e_k * F = t^{-C(k,2)} e_k(Y) F, in m=n variables.
#   sanity: nabla e_2 = s_2 + (q+t) s_11 (Haiman's convention).
# Negative controls must BREAK the claim.
import sys, itertools
from functools import reduce
import sympy as sp

q, t, s = sp.symbols('q t s')
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 4
P = sp.symbols('p1:%d' % (NMAX + 1))

def partitions(n, mx=None):
    if mx is None: mx = n
    if n == 0: yield (); return
    for k in range(min(n, mx), 0, -1):
        for r in partitions(n - k, k): yield (k,) + r

def conj(l): return tuple(sum(1 for x in l if x > j) for j in range(l[0])) if l else ()
def nn(l): return sum(i * x for i, x in enumerate(l))
def zee(l):
    from collections import Counter
    r = 1
    for k, m in Counter(l).items(): r *= k ** m * sp.factorial(m)
    return r
def pp(l): return reduce(lambda a, b: a * b, [P[i - 1] for i in l], sp.Integer(1))
def dominates(a, b):  # a >= b
    sa = sb = 0
    for i in range(max(len(a), len(b))):
        sa += a[i] if i < len(a) else 0; sb += b[i] if i < len(b) else 0
        if sa < sb: return False
    return True

def coeffs(F, n):  # p-basis coefficients
    F = sp.Poly(sp.expand(F), *P[:n])
    return {l: F.coeff_monomial(pp(l)) for l in partitions(n)}
def pair(F, G, n, w=lambda l: 1):
    a, b = coeffs(F, n), coeffs(G, n)
    return sp.together(sum(zee(l) * w(l) * a[l] * b[l] for l in partitions(n)))
def pleth(F, f):  # p_r -> f(r) p_r
    return sp.expand(F.subs({P[r - 1]: f(r) * P[r - 1] for r in range(1, NMAX + 1)}, simultaneous=True))

# monomial m_lam in p-basis: via h-basis duality is overkill; use n variables + linear solve
def m_in_p(n):
    xs = sp.symbols('x1:%d' % (n + 1))
    pv = {r: sum(x ** r for x in xs) for r in range(1, n + 1)}
    parts = list(partitions(n))
    def mono(l):
        lp = l + (0,) * (n - len(l))
        return sum(reduce(lambda a, b: a * b, [xs[i] ** e for i, e in enumerate(perm)], 1)
                   for perm in set(itertools.permutations(lp)))
    # matrix: p_mu = sum_lam M[mu][lam] m_lam ; coefficient of leading monomial of m_lam
    M = sp.zeros(len(parts))
    for i, mu in enumerate(parts):
        poly = sp.Poly(sp.expand(reduce(lambda a, b: a * b, [pv[r] for r in mu], 1)), *xs)
        for j, la in enumerate(parts):
            lp = la + (0,) * (n - len(la))
            M[i, j] = poly.coeff_monomial(reduce(lambda a, b: a * b, [xs[k] ** e for k, e in enumerate(lp)], 1))
    Minv = M.inv()  # m_lam = sum_mu Minv[lam][mu] p_mu
    return {la: sum(Minv[j, i] * pp(mu) for i, mu in enumerate(parts)) for j, la in enumerate(parts)}

def macdonaldP(n, qq, tt, mp):
    # Gram-Schmidt in monomial coordinates with a cancelled Gram matrix (no expression swell)
    parts = list(partitions(n))  # decreasing lex = linear ext of dominance; total for n<=5
    w = lambda l: reduce(lambda a, b: a * b, [(1 - qq ** k) / (1 - tt ** k) for k in l], 1)
    C = {la: coeffs(mp[la], n) for la in parts}           # m_la in p-coordinates
    G = {(a, b): sp.cancel(sum(zee(l) * w(l) * C[a][l] * C[b][l] for l in parts)) for a in parts for b in parts}
    V = {}                                                 # P_la as dict over m-basis
    def ip(u, v): return sp.cancel(sum(u[a] * v[b] * G[(a, b)] for a in u for b in v))
    for la in reversed(parts):
        u = {la: sp.Integer(1)}
        for mu in V:
            if dominates(la, mu) and mu != la:
                c = sp.cancel(ip({la: 1}, V[mu]) / ip(V[mu], V[mu]))
                for k, v in V[mu].items(): u[k] = sp.cancel(u.get(k, 0) - c * v)
        V[la] = u
    return {la: sp.expand(sum(c * mp[k] for k, c in V[la].items())) for la in parts}

def cnu(l, qq, tt):
    lc = conj(l); r = 1
    for i, row in enumerate(l):
        for j in range(row):
            a = row - j - 1; lg = lc[j] - i - 1
            r *= (1 - qq ** a * tt ** (lg + 1))
    return r

def schur(l, n):  # Jacobi-Trudi with h_k = sum_{|mu|=k} p_mu/z_mu
    h = lambda k: 0 if k < 0 else (1 if k == 0 else sum(pp(mu) / zee(mu) for mu in partitions(k)))
    L = len(l)
    return sp.expand(sp.Matrix(L, L, lambda i, j: h(l[i] - i + j)).det())

def simp(e): return sp.simplify(sp.factor(sp.together(e)))
def iszero(F, n):
    return all(simp(c) == 0 for c in coeffs(F, n).values())

def gh_Htilde(mu, n, S):
    # unknown combination of Schur functions, solve GH triangularities + <H, s_n> = 1
    parts = list(partitions(n))
    cs = sp.symbols('c0:%d' % len(parts))
    H = sum(c * S[l] for c, l in zip(cs, parts))
    eqs = []
    A = pleth(H, lambda r: 1 - q ** r); B = pleth(H, lambda r: 1 - t ** r)
    for l in parts:
        if not dominates(l, mu): eqs.append(pair(A, S[l], n))
        if not dominates(l, conj(mu)): eqs.append(pair(B, S[l], n))
    eqs.append(pair(H, S[(n,)], n) - 1)
    sol = sp.solve([sp.numer(sp.together(e)) for e in eqs], cs, dict=True)
    assert len(sol) == 1 and all(c in sol[0] for c in cs), ('GH not unique', mu)
    return sp.expand(H.subs(sol[0]))

def expand_in(F, basis, n):  # solve F = sum a_l basis[l]
    parts = list(partitions(n))
    a = sp.symbols('a0:%d' % len(parts))
    G = F - sum(ai * basis[l] for ai, l in zip(a, parts))
    sol = sp.solve([sp.together(c) for c in coeffs(G, n).values()], a, dict=True)[0]
    return {l: sp.factor(sol[ai]) for ai, l in zip(a, parts)}

eps_sign = lambda r: (-1) ** (r - 1)          # p_r[-eps X] = (-1)^{r-1} p_r   (printed)
phi = lambda F, tt=t: pleth(F, lambda r: eps_sign(r) / (1 - tt ** r))
phiinv = lambda F, tt=t: pleth(F, lambda r: eps_sign(r) * (1 - tt ** r))

# ---------------- printed Sec. 2 operators (for P4) -----------------
def star_e(lam, m):
    xs = sp.symbols('y1:%d' % (m + 1))
    def sw(f, i): return f.subs({xs[i]: xs[i + 1], xs[i + 1]: xs[i]}, simultaneous=True)
    def T(f, i):  # 0-based i <-> T_{i+1}
        g = sw(f, i)
        return sp.expand(t * g + sp.cancel((t - 1) * xs[i + 1] * (g - f) / (xs[i] - xs[i + 1])))
    def Tinv(f, i): return sp.expand(T(f, i) / t - (1 - 1 / t) * f)  # (R4)/(H): T^-1 = t^-1 T - (1 - t^-1)
    def pi(f): return sp.expand(xs[0] * f.subs({**{xs[j]: xs[j + 1] for j in range(m - 1)}, xs[m - 1]: s * xs[0]}, simultaneous=True))
    def Y(f, i):  # 1-based i: t^{m-i} T_{i-1}..T_1 pi T_{m-1}^{-1}..T_i^{-1}, rightmost first
        for j in range(i, m): f = Tinv(f, j - 1)  # T_i^{-1} acts first
        f = pi(f)
        for j in range(1, i): f = T(f, j - 1)    # then T_1, ..., T_{i-1}
        return sp.expand(t ** (m - i) * f)
    def ekY(F, k):
        out = 0
        for B in itertools.combinations(range(1, m + 1), k):
            g = F
            for b in reversed(B): g = Y(g, b)   # product Y_{b1}...Y_{bk}, left to right
            out += g
        return sp.expand(out)
    F = sp.Integer(1)
    for k in reversed(lam): F = sp.expand(ekY(F, k) / t ** sp.binomial(k, 2))
    return F, xs

def to_poly_in_x(F, xs):
    return sp.expand(F.subs({P[r - 1]: sum(x ** r for x in xs) for r in range(1, NMAX + 1)}, simultaneous=True))

# ---------------- run -----------------
results = []
def rep(name, ok):
    results.append((name, ok)); print(('PASS ' if ok else 'FAIL ') + name, flush=True)

for n in range(1, NMAX + 1):
    parts = list(partitions(n))
    mp = m_in_p(n)
    S = {l: schur(l, n) for l in parts}
    Pq = macdonaldP(n, q, t, mp)              # P(x;q,t)
    Jq = {l: sp.expand(cnu(l, q, t) * Pq[l]) for l in parts}
    Jqi = {l: sp.expand(sp.together(Jq[l].subs(t, 1 / t))) for l in parts}  # J(x;q,1/t)
    Ht = {l: gh_Htilde(l, n, S) for l in parts}
    for l in parts:
        # (P1) printed normalisation of Ht vs GH
        rhs = t ** nn(l) * pleth(Jqi[l], lambda r: 1 / (1 - t ** (-r)))
        rep('P1 n=%d nu=%s Ht = t^n(nu) J[X/(1-1/t);q,1/t]' % (n, l), iszero(Ht[l] - rhs, n))
        # (P2) phi^{-1} Ht(x;s,t) = (-1)^n t^{n+n(nu)} J(x;s,1/t)
        lhs = phiinv(Ht[l].subs(q, s))
        J_s = Jqi[l].subs(q, s)
        rep('P2 n=%d nu=%s phi^-1 Ht = (-1)^n t^(n+n(nu)) J(s,1/t)' % (n, l), iszero(lhs - (-1) ** n * t ** (n + nn(l)) * J_s, n))
        # controls for P2 (must FAIL)
        rep('  ctrl-sign(+1 for odd n) breaks' if n % 2 else '  ctrl-sign(-(..)) breaks', not iszero(lhs + (-1) ** n * t ** (n + nn(l)) * J_s, n) or False)
        if conj(l) != l:
            rep('  ctrl n(nu)->n(nu\') breaks', not iszero(lhs - (-1) ** n * t ** (n + nn(conj(l))) * J_s, n))
        if n >= 2:  # phi without -eps: phi^{-1}(G)=G[X(1-t)]; result must NOT be a multiple of J(s,1/t)
            lhs_noeps = pleth(Ht[l].subs(q, s), lambda r: (1 - t ** r))
            ratios = set(simp(c1 / c2) for c1, c2 in zip(coeffs(lhs_noeps, n).values(), coeffs(J_s, n).values()) if c2 != 0)
            zeros_ok = all((c1 == 0) == (c2 == 0) for c1, c2 in zip(coeffs(lhs_noeps, n).values(), coeffs(J_s, n).values()))
            rep('  ctrl phi without -eps: not prop. to J (breaks)', len(ratios) > 1 or not zeros_ok)
    if n == 2:
        Hs = Ht
        ev = {l: t ** nn(l) * q ** nn(conj(l)) for l in parts}
        e2 = pp((1, 1)) / 2 - pp((2,)) / 2
        a = expand_in(e2, Ht, n)
        nab = sp.expand(sum(a[l] * ev[l] * Ht[l] for l in parts))
        rep('sanity nabla e_2 = s_2 + (q+t) s_11', iszero(nab - S[(2,)] - (q + t) * S[(1, 1)], n))
    # (P3) N = phi^{-1} nabla_{q=s} phi on e_lam basis
    Ps = {l: Pq[l].subs({q: s}).subs(t, 1 / t) for l in parts}   # P(x;s,1/t)
    Hs = {l: Ht[l].subs(q, s) for l in parts}
    eN = {l: t ** nn(l) * s ** nn(conj(l)) for l in parts}
    def ebasis(l):
        e = lambda k: sum((-1) ** (k - len(mu)) * pp(mu) / zee(mu) for mu in partitions(k))
        return sp.expand(reduce(lambda a, b: a * b, [e(k) for k in l], 1))
    def Nop(F, evf=eN):
        a = expand_in(F, Ps, n); return sp.expand(sum(a[l] * evf[l] * Ps[l] for l in parts))
    def conj_nabla(F, evf, phif=phi, phiif=phiinv, H=Hs):
        a = expand_in(phif(F), H, n); return phiif(sp.expand(sum(a[l] * evf[l] * H[l] for l in parts)))
    for lam in parts:
        E = ebasis(lam)
        NE = Nop(E)
        rep('P3 n=%d lam=%s N e = phi^-1 nabla phi e' % (n, lam), iszero(NE - conj_nabla(E, eN), n))
        if lam == parts[0] or lam == parts[-1]:
            continue  # e_n and e_1^n handled below via generic controls
    # controls for P3 on the whole space: wrong eigenvalue, s<->t, plethystic sign
    bad_ev = {l: t ** nn(conj(l)) * s ** nn(l) for l in parts}
    swap_ev = {l: s ** nn(l) * t ** nn(conj(l)) for l in parts}
    Hswap = {l: Ht[l].subs({q: t, t: s}, simultaneous=True) for l in parts}  # nabla at (q,t)=(t,s)
    ev_swap = {l: s ** nn(l) * t ** nn(conj(l)) for l in parts}
    phi_plus = lambda F, tt=t: pleth(F, lambda r: 1 / (1 - tt ** r))
    phii_plus = lambda F, tt=t: pleth(F, lambda r: (1 - tt ** r))
    if n >= 2:
        def breaks(f): return any(not iszero(Nop(ebasis(l)) - f(ebasis(l)), n) for l in parts)
        rep('  ctrl n(nu)<->n(nu\') eigenvalue breaks (n=%d)' % n, breaks(lambda F: conj_nabla(F, bad_ev)))
        sig = lambda G: G.subs({s: t, t: s}, simultaneous=True)
        rep('  ctrl s<->t (whole RHS operator) breaks (n=%d)' % n, breaks(lambda F: sig(conj_nabla(F, eN))))
        rep('  ctrl phi without -eps breaks (n=%d)' % n, breaks(lambda F: conj_nabla(F, eN, phi_plus, phii_plus)))
    # (P4) end-to-end with printed Sec.2 operators, m = n variables
    if n <= NMAX:
        m = n
        for lam in parts:
            if lam[0] > m: continue
            Fst, xs = star_e(lam, m)
            rhs = to_poly_in_x(conj_nabla(ebasis(lam), eN), xs)
            rep('P4 n=%d lam=%s t^n(lam\') e*_lam = phi^-1 nabla phi e_lam (m=%d, printed T,pi,Y)' % (n, lam, m),
                sp.simplify(t ** nn(conj(lam)) * Fst - rhs) == 0)
            if lam != conj(lam):
                rep('  ctrl P4 with t^n(lam) instead of t^n(lam\') breaks', sp.simplify(t ** nn(lam) * Fst - rhs) != 0)

print('ALL', all(ok for _, ok in results), len(results))
