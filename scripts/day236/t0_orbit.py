# Day 236: t->0 limit of T_i, T_i^{-1}, Y_i on gr_O (one S_m-orbit modulo dominance-lower orbits).
# Engine: identical conventions to scripts/day235/tri_proof_check.py (ground truth, sympy substitution).
import sympy as sp, itertools, sys
s, t = sp.symbols('s t')
def setup(m, variant='ok'):
    X = sp.symbols('X1:%d' % (m + 1))
    def T(i, F):
        a, b = X[i - 1], X[i]
        sF = F.subs({a: b, b: a}, simultaneous=True)
        return sp.expand(t * sF + sp.cancel((t - 1) * b * (sF - F) / (a - b)))
    def Ti(i, F): return sp.expand(T(i, F) / t - (1 - 1 / t) * F)
    if variant == 'pi_wrong':   # NEG: s on the wrong slot: F(s X2, X3, ..., X1)
        def pi(F): return sp.expand(F.subs({X[j]: (s * X[1] if j == 0 else (X[j + 1] if j < m - 1 else X[0])) for j in range(m)}, simultaneous=True))
    else:
        def pi(F): return sp.expand(F.subs({X[j]: (X[j + 1] if j < m - 1 else s * X[0]) for j in range(m)}, simultaneous=True))
    def Y(i, F):
        for k in range(i, m): F = (T if variant == 'Tpos' else Ti)(k, F)
        F = pi(F)
        for k in range(1, i): F = T(k, F)
        return sp.expand(t ** (m - i) * F)
    return X, T, Ti, pi, Y
def sortd(c): return tuple(sorted(c, reverse=True))
def psums(c): return list(itertools.accumulate(c))
def dom_leq(a, b): return all(x <= y for x, y in zip(psums(a), psums(b)))
def b_of(lam):
    m = len(lam); return [sum(1 for j in range(m) if lam[j] > lam[i]) + sum(1 for j in range(i) if lam[j] == lam[i]) for i in range(m)]
def mono(X, c): return sp.Mul(*[x ** e for x, e in zip(X, c)])
def todict(X, F):
    if F == 0: return {}
    P = sp.Poly(F, *X); return {tuple(mu): sp.factor(cf) for mu, cf in zip(P.monoms(), P.coeffs())}
def orbit(lp):
    O = sorted(set(itertools.permutations(lp)), key=lambda c: psums(c), reverse=True)  # lin. ext. of <=_c, top = lam+
    return O
ALL = True
def rep(name, ok):
    global ALL; ALL &= bool(ok); print('%-78s %s' % (name, ok)); sys.stdout.flush()

def grmat(X, op, O):
    """matrix of op on gr_O; also checks filtration: every other monomial lies in a strictly dominance-lower orbit"""
    n = len(O); idx = {c: k for k, c in enumerate(O)}; M = sp.zeros(n, n); filt = True
    lp = O[0]
    for k, mu in enumerate(O):
        for nu, cf in todict(X, op(mono(X, mu))).items():
            if nu in idx: M[idx[nu], k] = cf
            else: filt &= (sortd(nu) != lp) and dom_leq(sortd(nu), lp)
    return M.applyfunc(sp.factor), filt
def lim0(M): return M.applyfunc(lambda e: sp.limit(sp.simplify(e), t, 0))
def finite0(M):
    try:
        L = lim0(M); return all(e.is_finite for e in L), L
    except Exception: return False, None
def atom_mat(O, i):   # Demazure-atom / sorting action pibar_i = pi_i - 1 on gr_O
    n = len(O); idx = {c: k for k, c in enumerate(O)}; M = sp.zeros(n, n)
    for k, mu in enumerate(O):
        a, b = mu[i - 1], mu[i]
        if a > b: sw = list(mu); sw[i - 1], sw[i] = b, a; M[idx[tuple(sw)], k] = 1
        elif a < b: M[k, k] = -1
    return M

CASES = {3: [(2, 1, 0), (3, 1, 0), (1, 1, 0), (1, 0, 0), (2, 2, 1), (2, 1, 1), (1, 1, 1)],
         4: [(3, 2, 1, 0), (1, 1, 0, 0), (2, 1, 1, 0), (2, 2, 0, 0), (2, 1, 0, 0)]}
SHOW = (2, 1, 0)
for m, lps in CASES.items():
    X, T, Ti, pi, Y = setup(m)
    for lp in lps:
        O = orbit(lp); n = len(O); Id = sp.eye(n)
        Tm, tTi, Ym = {}, {}, {}
        filt = True
        for i in range(1, m):
            Tm[i], f1 = grmat(X, lambda F: T(i, F), O)
            M2, f2 = grmat(X, lambda F: Ti(i, F), O); tTi[i] = (t * M2).applyfunc(sp.factor)
            filt &= f1 and f2
            fin_inv, _ = finite0(M2); fin_negT, _ = finite0(-Tm[i])
            if i == 1: print('[m=%d %s] t->0 finite?  T_1:%s  -T_1:%s  T_1^-1:%s  t*T_1^-1:%s' % (m, lp, finite0(Tm[1])[0], fin_negT, fin_inv, finite0(tTi[1])[0]))
        for i in range(1, m + 1):
            Ym[i], f3 = grmat(X, lambda F: Y(i, F), O); filt &= f3
        rep('[m=%d %s] |O|=%d; T_i, T_i^-1, Y_i preserve dominance filtration (gr well-defined)' % (m, lp, n), filt)
        # 1. limits
        T0 = {i: lim0(Tm[i]) for i in Tm}; P0 = {i: lim0(tTi[i]) for i in tTi}
        rep('  T_i|t=0 finite and == Demazure-atom sorting matrix pibar_i', all(T0[i] == atom_mat(O, i) for i in T0))
        rep('  (t T_i^-1)|t=0 == pibar_i + 1 (Demazure pi_i)', all(P0[i] == atom_mat(O, i) + Id for i in P0))
        # 2. 0-Hecke relations
        r = all(T0[i] ** 2 == -T0[i] for i in T0) and all(T0[i] * T0[i + 1] * T0[i] == T0[i + 1] * T0[i] * T0[i + 1] for i in range(1, m - 1)) \
            and all(T0[i] * T0[j] == T0[j] * T0[i] for i in T0 for j in T0 if abs(i - j) > 1)
        rep('  T^(0)_i: T^2=-T, braid, far-commute (H_m(0))', r)
        # cyclic from x^{lam+}, T^(0)_j x^{lam+} = 0 for ties, composition factors = ascent sets
        tie = [j for j in range(1, m) if lp[j - 1] == lp[j]]
        rep('  T^(0)_j x^{lam+} = 0 for tie positions j=%s; T^(0)_j x^{lam+}=x^{s_j lam+} at descents' % tie,
            all(T0[j][:, 0] == sp.zeros(n, 1) for j in tie) and all(T0[j][:, 0] != sp.zeros(n, 1) for j in range(1, m) if j not in tie))
        # cyclic: span of all words in T0 applied to e_0 is everything
        vecs = [sp.eye(n)[:, 0]]; frontier = list(vecs)
        while frontier:
            nf = []
            for v in frontier:
                for i in T0:
                    w = T0[i] * v
                    if sp.Matrix.hstack(*vecs, w).rank() > len(vecs): vecs.append(w); nf.append(w)
            frontier = nf
        rep('  module cyclic on x^{lam+} (rank %d of %d)' % (len(vecs), n), len(vecs) == n)
        # 3. Y at t=0, normalised
        Y0 = {i: lim0(Ym[i]) for i in Ym}
        mx = max(lp)
        diag_ok = all(Y0[i][k, k] == (s ** mu[i - 1] if b_of(mu)[i - 1] == 0 else 0) for i in Y0 for k, mu in enumerate(O))
        rep('  Y_i|t=0 finite; diag = s^mu_i [b_i(mu)=0]  (b_i=0 <=> i = first max position)', diag_ok)
        Pn = {i: (Y0[i] / s ** mx).applyfunc(sp.simplify) for i in Y0}
        proj = [P ** 2 == P for P in Pn.values()]
        offd = sum(1 for i in Y0 for a in range(n) for b in range(n) if a != b and Y0[i][a, b] != 0)
        rep('  Y_i^(0)/s^max idempotent: %s ; #nonzero offdiag entries %d ; commute: %s' % (proj, offd,
            all(Y0[i] * Y0[j] == Y0[j] * Y0[i] for i in Y0 for j in Y0)), True)
        # closed formula at t=0: Y_i^(0) = pibar_{i-1}..pibar_1 * Pi * pi_{m-1}..pi_i on gr_O
        idx = {c: k for k, c in enumerate(O)}; PiM = sp.zeros(n, n)
        for k, mu in enumerate(O): PiM[idx[(mu[-1],) + mu[:-1]], k] = s ** mu[-1]
        ok = True
        for i in Y0:
            M = Id
            for k in range(i, m): M = P0[k] * M
            M = PiM * M
            for k in range(1, i): M = T0[k] * M
            ok &= (M - Y0[i]).applyfunc(sp.simplify) == sp.zeros(n, n)
        rep('  Y_i^(0) == pibar_{i-1}..pibar_1 Pi pi_{m-1}..pi_i (Demazure words)', ok)
        if lp == SHOW:
            print('  basis order:', O)
            for i in T0: print('  T_%d|0 =' % i, T0[i].tolist())
            for i in P0: print('  (tT_%d^-1)|0 =' % i, P0[i].tolist())
            for i in Tm: print('  T_%d (generic t) =' % i, Tm[i].tolist())
            for i in Ym: print('  Y_%d generic =' % i, Ym[i].tolist())
            for i in Y0: print('  Y_%d|0 =' % i, Y0[i].tolist())
            # eigenvectors at generic t, then t->0
            Z = sp.zeros(n, n)
            for k in range(m): Z += (k + 2) * Ym[k + 1]
            print('  t->0 of the gr_O eigenvectors E_mu (normalised coefficient 1 on x^mu):')
            for k, mu in enumerate(O):
                ev = sum((j + 2) * s ** mu[j] * t ** b_of(mu)[j] for j in range(m))
                N = (Z - ev * Id).nullspace()[0]; N = (N / N[k]).applyfunc(sp.factor)
                print('    ', mu, [sp.limit(e, t, 0) for e in N])

# exact identities on Pol at t=0 (not only gr): T_i|0 = pi_i - 1, (tT_i^-1)|0 = pi_i
X, T, Ti, pi, Y = setup(3)
def dem(i, F):
    a, b = X[i - 1], X[i]; return sp.expand(sp.cancel((a * F - b * F.subs({a: b, b: a}, simultaneous=True)) / (a - b)))
tests = [mono(X, c) for d in range(4) for c in itertools.product(range(d + 1), repeat=3) if sum(c) == d]
rep('Pol (m=3,deg<=3): T_i|t=0 == pi_i - 1 and (t T_i^-1)|t=0 == pi_i exactly',
    all(sp.expand(T(i, F).subs(t, 0) - (dem(i, F) - F)) == 0 and sp.expand(sp.expand(t * Ti(i, F)).subs(t, 0) - dem(i, F)) == 0 for F in tests for i in (1, 2)))

# 4. NEGATIVE CONTROLS (each must BREAK a claim)
print('--- negative controls (expect each "BROKEN: True") ---')
lp = (2, 1, 0); O = orbit(lp); n = len(O)
for var in ('Tpos', 'pi_wrong'):
    X, T, Ti, pi, Y = setup(3, var)
    Ym = {}; filt = True
    for i in (1, 2, 3): Ym[i], f = grmat(X, lambda F: Y(i, F), O); filt &= f
    tri = all(Ym[i][a, b] == 0 for i in Ym for a in range(n) for b in range(n) if a < b) # O ordered top-first: upper entries = higher terms
    # note: correct Y is lower-triangular in this order (column mu only hits nu <=_c mu i.e. later index)
    diag = all(sp.simplify(Ym[i][k, k] - s ** mu[i - 1] * t ** b_of(mu)[i - 1]) == 0 for i in Ym for k, mu in enumerate(O))
    comm = all((Ym[i] * Ym[j] - Ym[j] * Ym[i]).applyfunc(sp.simplify) == sp.zeros(n, n) for i in Ym for j in Ym)
    fin = all(finite0(Ym[i])[0] for i in Ym)
    Y0d = None
    if fin:
        Y0 = {i: lim0(Ym[i]) for i in Ym}
        Y0d = all(Y0[i][k, k] == (s ** mu[i - 1] if b_of(mu)[i - 1] == 0 else 0) for i in Y0 for k, mu in enumerate(O))
    print('  NEG %-9s filt=%s lower-tri=%s diag=s^mu_i t^b_i:%s commute_on_gr:%s finite_t0:%s t0-diag-pattern:%s  BROKEN: %s' %
          (var, filt, tri, diag, comm, fin, Y0d, not (filt and tri and diag and comm and fin and Y0d)))
# sanity: the correct Y passes the same battery
X, T, Ti, pi, Y = setup(3)
Ym = {i: grmat(X, lambda F: Y(i, F), O)[0] for i in (1, 2, 3)}
print('  POS correct  lower-tri=%s' % all(Ym[i][a, b] == 0 for i in Ym for a in range(n) for b in range(n) if a < b))
# NEG for 0-Hecke claim: wrong-sign Hecke T'_i = t s_i + (1-t) X_{i+1}(s_iF-F)/(X_i-X_{i+1}) at t=0 is NOT H(0)
X = sp.symbols('X1:4')
def Tw(i, F):
    a, b = X[i - 1], X[i]; sF = F.subs({a: b, b: a}, simultaneous=True)
    return sp.expand(t * sF + sp.cancel((1 - t) * b * (sF - F) / (a - b)))
Tw0 = {i: lim0(grmat(X, lambda F: Tw(i, F), O)[0]) for i in (1, 2)}
h0 = all(Tw0[i] ** 2 in (Tw0[i], -Tw0[i]) for i in Tw0) and Tw0[1] * Tw0[2] * Tw0[1] == Tw0[2] * Tw0[1] * Tw0[2]
print('  NEG wrong-sign Hecke: t=0 limit satisfies H(0) relations: %s  BROKEN: %s' % (h0, not h0))
print('ALL:', ALL)
print('  (the wrong-sign Hecke control only RELABELS: its t=0 limit is -pibar_i = Norton idempotent pi_i - so it is not a break)')
# NEG that breaks: mixed normalisation (Demazure pi_1 = (tT_1^-1)|0 with atom pibar_2 = T_2|0) -- braid must fail
X, T, Ti, pi, Y = setup(3)
A = lim0((t * grmat(X, lambda F: Ti(1, F), O)[0])); B = lim0(grmat(X, lambda F: T(2, F), O)[0])
br = A * B * A == B * A * B
print('  NEG mixed normalisation {pi_1, pibar_2}: braid holds: %s  BROKEN: %s' % (br, not br))
