# Day 235: cross-checks for the hand proofs of (A) presentation and (B) triangularity.
# Ground-truth operators: sympy substitution engine, identical conventions to scripts/day233/tri_check.py.
import sympy as sp, itertools, sys
s, t = sp.symbols('s t')
def setup(m, variant='ok'):
    X = sp.symbols('X1:%d' % (m + 1))
    def T(i, F):
        a, b = X[i - 1], X[i]
        sF = F.subs({a: b, b: a}, simultaneous=True)
        return sp.expand(t * sF + sp.cancel((t - 1) * b * (sF - F) / (a - b)))
    def Ti(i, F): return sp.expand(T(i, F) / t - (1 - 1 / t) * F)
    def pi(F): return sp.expand(F.subs({X[j]: (X[j + 1] if j < m - 1 else s * X[0]) for j in range(m)}, simultaneous=True))
    def Y(i, F):
        for k in range(i, m): F = (Ti if variant == 'ok' else T)(k, F)   # NEG 'Tpos': T instead of T^{-1}
        F = pi(F)
        for k in range(1, i): F = T(k, F)
        return sp.expand(t ** (m - i) * F)
    return X, T, Ti, pi, Y
def comps(m, d): return [c for c in itertools.product(range(d + 1), repeat=m) if sum(c) == d]
def psums(c): return list(itertools.accumulate(c))
def dom_leq(a, b): return all(x <= y for x, y in zip(psums(a), psums(b)))
def sortd(c): return tuple(sorted(c, reverse=True))
def b_of(lam):
    m = len(lam); return [sum(1 for j in range(m) if lam[j] > lam[i]) + sum(1 for j in range(i) if lam[j] == lam[i]) for i in range(m)]
def mono(X, c): return sp.Mul(*[x ** e for x, e in zip(X, c)])
def todict(X, F):
    P = sp.Poly(F, *X); return {tuple(mu): sp.factor(cf) for mu, cf in zip(P.monoms(), P.coeffs())}
ALL = True
def rep(name, ok):
    global ALL; ALL &= ok; print('%-70s %s' % (name, ok)); sys.stdout.flush()

# 1. T on monomials (m=2) vs hand formula
X, T, Ti, pi, Y = setup(2)
ok = True
for a in range(6):
    for b in range(6):
        F = mono(X, (a, b)); got = T(1, F); gti = Ti(1, F)
        btw = sum(mono(X, (c, a + b - c)) for c in range(min(a, b) + 1, max(a, b)))
        if a > b: pred = mono(X, (b, a)) + (1 - t) * btw; predi = mono(X, (b, a)) / t - (1 - 1 / t) * F + (1 - t) / t * btw
        elif a < b: pred = t * mono(X, (b, a)) + (t - 1) * F + (t - 1) * btw; predi = mono(X, (b, a)) + (1 - 1 / t) * btw
        else: pred = t * F; predi = F / t
        ok &= sp.simplify(got - pred) == 0 and sp.simplify(gti - predi) == 0
rep('1. T, T^-1 on x^(a,b), a,b<=5, match hand formulas', ok)

# 2. braid relation on the 6 Artin monomials X2^a X3^b (a<=1,b<=2), m=3; quadratic on 1, X2 (m=2)
X, T, Ti, pi, Y = setup(3)
ok = all(sp.expand(T(1, T(2, T(1, mono(X, (0, a, b))))) - T(2, T(1, T(2, mono(X, (0, a, b)))))) == 0 for a in range(2) for b in range(3))
rep('2a. braid T1T2T1=T2T1T2 on Artin basis X2^a X3^b (a<=1,b<=2)', ok)
ok = all(sp.expand(T(1, T(1, F)) - (t - 1) * T(1, F) - t * F) == 0 for F in (sp.Integer(1), X[1], X[0], X[0]**2*X[1]))
rep('2b. quadratic (T-t)(T+1)=0 on 1, X2 (+ X1, X1^2X2 extra)', ok)
# spanning claim for the Artin basis: X1, X2^2, X3^3 expressed via e's (sanity, deg 3 monomials)
e1, e2, e3 = X[0] + X[1] + X[2], X[0]*X[1] + X[0]*X[2] + X[1]*X[2], X[0]*X[1]*X[2]
ok = sp.expand(X[2]**3 - (e1*X[2]**2 - e2*X[2] + e3)) == 0 and sp.expand(X[1]**2 - ((e1 - X[2])*X[1] - (e2 - X[2]*(e1 - X[2])))) == 0
rep('2c. Artin reduction identities x3^3, x2^2 (hand proof of spanning)', ok)

# 3. (c),(d) and C1 at m=5 small degree (new range; m=3,4 done Day 233)
for m, D in ((5, 2),):
    X, T, Ti, pi, Y = setup(m)
    tests = [mono(X, c) for d in range(D + 1) for c in comps(m, d)]
    ok = all(sp.expand(pi(T(i, F)) - T(i + 1, pi(F))) == 0 for F in tests for i in range(1, m - 1))
    rep('3a. m=%d pi T_i = T_{i+1} pi, deg<=%d' % (m, D), ok)
    ok = all(sp.expand(pi(pi(T(m - 1, F))) - T(1, pi(pi(F)))) == 0 for F in tests)
    rep('3b. m=%d pi^2 T_{m-1} = T_1 pi^2, deg<=%d' % (m, D), ok)
    ok = all(sp.expand(Y(i, Y(j, F)) - Y(j, Y(i, F))) == 0 for F in tests for i in range(1, m + 1) for j in range(i + 1, m + 1))
    rep('3c. m=%d C1 Y_iY_j=Y_jY_i, deg<=%d' % (m, D), ok)

# 4. triangularity in the PROOF's order: mu+ < lam+ strictly, or same orbit with partial sums mu <= lam (mu != lam);
#    diagonal s^lam_i t^b_i(lam); and the gr-level path model reproduces the same-orbit part EXACTLY.
def path_model(lam, i):
    m = len(lam); out = {}
    # phase 1: T_j^{-1}, j=i..m-1 acting on gr (mod lower dominance)
    states = {tuple(lam): sp.Integer(1)}
    def apply_gr(states, j, inv):
        new = {}
        for c, cf in states.items():
            p, q = c[j - 1], c[j]; sw = list(c); sw[j - 1], sw[j] = q, p; sw = tuple(sw)
            if inv:
                terms = [(c, 1 / t)] if p == q else ([(sw, 1)] if p < q else [(sw, 1 / t), (c, -(1 - 1 / t))])
            else:
                terms = [(c, t)] if p == q else ([(sw, 1)] if p > q else [(sw, t), (c, t - 1)])
            for d, k in terms: new[d] = new.get(d, 0) + cf * k
        return new
    for j in range(i, m): states = apply_gr(states, j, True)
    states = {(c[-1],) + c[:-1]: cf * s ** c[-1] for c, cf in states.items()}
    for j in range(1, i): states = apply_gr(states, j, False)
    return {c: sp.factor(cf * t ** (m - i)) for c, cf in states.items() if sp.simplify(cf) != 0}
for m, D in ((2, 5), (3, 4), (4, 3), (5, 2)):
    X, T, Ti, pi, Y = setup(m)
    okorder = okdiag = okpath = True; nsame = nlow = 0; inv_ok = True
    for d in range(D + 1):
        for lam in comps(m, d):
            for i in range(1, m + 1):
                Dd = todict(X, Y(i, mono(X, lam)))
                okdiag &= sp.simplify(Dd.get(lam, 0) - s ** lam[i - 1] * t ** b_of(lam)[i - 1]) == 0
                same = {}
                for mu, cf in Dd.items():
                    if mu == lam: continue
                    if sortd(mu) == sortd(lam):
                        nsame += 1; same[mu] = cf
                        okorder &= dom_leq(mu, lam)
                    else:
                        nlow += 1; okorder &= dom_leq(sortd(mu), sortd(lam))
                same[lam] = Dd.get(lam, 0)
                pm = path_model(lam, i)
                keys = set(pm) | set(same)
                okpath &= all(sp.simplify(pm.get(k, 0) - same.get(k, 0)) == 0 for k in keys)
    rep('4. m=%d d<=%d: diag = s^lam_i t^b_i' % (m, D), okdiag)
    rep('4. m=%d d<=%d: off-diag lower in proof order (%d same-orbit, %d lower-dom)' % (m, D, nsame, nlow), okorder)
    rep('4. m=%d d<=%d: gr path model == actual same-orbit part' % (m, D), okpath)

# 5. NEGATIVE CONTROLS (must FAIL)
# 5a. reversed within-orbit order (mu >= lam in partial sums) must be violated
X, T, Ti, pi, Y = setup(3)
viol = False
for lam in comps(3, 3):
    for i in (1, 2, 3):
        for mu in todict(X, Y(i, mono(X, lam))):
            if mu != lam and sortd(mu) == sortd(lam) and not dom_leq(lam, mu): viol = True
print('5a. NEG reversed orbit order violated (expect True):', viol)
# 5b. Y with T instead of T^{-1} in phase 1: commutativity must fail
X, T, Ti, pi, Yb = setup(3, 'Tpos')
F = mono(X, (1, 0, 1))
fail = any(sp.expand(Yb(i, Yb(j, F)) - Yb(j, Yb(i, F))) != 0 for i in (1, 2, 3) for j in (1, 2, 3))
print('5b. NEG T-instead-of-T^-1 Y commutativity fails (expect True):', fail)
# 5c. the PROVE.md misstatement pi x^lam = s^lam_1 x^(lam_2..lam_m,lam_1) must disagree with the engine
X, T, Ti, pi, Y = setup(3)
wrong = sp.expand(s**2 * mono(X, (0, 1, 2))); right = pi(mono(X, (2, 0, 1)))
print('5c. NEG PROVE.md pi formula disagrees with engine (expect True):', sp.expand(wrong - right) != 0, ' engine pi x^(2,0,1) =', right)
ALL &= viol and fail
print('ALL:', ALL)
