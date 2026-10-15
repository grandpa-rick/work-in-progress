"""Day 240: check Clio's Y-match (review 2026-10-14, item 1 / F3).
Operators implemented from the PRINTED definitions:
 DFK (arXiv 1704.00154, deformDAHA.tex sec 'functional representation', l.~127-133; Def. of Y_{i,n} l.~180):
   rho(T_i) = theta s_i + (theta-theta^-1)/(x_i x_{i+1}^-1 - 1) (s_i - 1);  rho(T_i^-1)=rho(T_i)-theta+theta^-1
   rho(pi) f = f(x_2,...,x_N, q^-1 x_1);  Y_i = T_i...T_{N-1} pi^-1 T_1^-1...T_{i-1}^-1
   Y_{i,n} = (X_1..X_{i-1})^-n Y_i (X_1..X_i)^n
 Hikita (arXiv 2503.23597, eqs (1)-(3), p.10-11):
   T_i F = t s_i F + (t-1)(F - s_i F)/(1 - X_i X_{i+1}^-1)
   T_i^-1 F = s_i F + (1-t^-1) (X_i X_{i+1}^-1 F - s_i F)/(1 - X_i X_{i+1}^-1)   (proof of Lemma 2.11)
   Pi F = X_1 F(X_2,...,X_m, q^-1 X_1)
   Y_i = t^{m-i} T_{i-1}...T_1 Pi T_{m-1}^-1 ... T_i^-1
Words act as operators: rightmost factor applied first.
"""
import sympy as sp, itertools, sys
th, q = sp.symbols('theta q')

def mk(N):
    return sp.symbols('x1:%d' % (N+1))

def simp(e): return sp.factor(sp.cancel(sp.together(e)))

def swap(f, x, i):  # s_i, i 1-based
    a, b = x[i-1], x[i]
    return f.subs({a: b, b: a}, simultaneous=True)

def sub(f, x, new): return f.subs(dict(zip(x, new)), simultaneous=True)

# ---------- DFK ----------
def D_T(f, x, i, th=th):
    s = swap(f, x, i)
    return sp.cancel(th*s + (th-1/th)/(x[i-1]/x[i]-1)*(s-f))
def D_Ti(f, x, i, th=th):
    return sp.cancel(D_T(f, x, i, th) - th*f + f/th)
def D_pi(f, x):
    N = len(x); return sub(f, x, list(x[1:]) + [x[0]/q])
def D_pii(f, x):
    N = len(x); return sub(f, x, [q*x[-1]] + list(x[:-1]))
def D_Y(f, x, i, th=th):   # T_i..T_{N-1} pi^-1 T_1^-1..T_{i-1}^-1
    N = len(x)
    for j in range(i-1, 0, -1): f = D_Ti(f, x, j, th)
    f = D_pii(f, x)
    for j in range(N-1, i-1, -1): f = D_T(f, x, j, th)
    return sp.cancel(f)
def D_Yinv(f, x, i, th=th):  # T_{i-1}..T_1 pi T_{N-1}^-1..T_i^-1
    N = len(x)
    for j in range(i, N): f = D_Ti(f, x, j, th)
    f = D_pi(f, x)
    for j in range(1, i): f = D_T(f, x, j, th)
    return sp.cancel(f)
def prod(xs):
    r = sp.Integer(1)
    for v in xs: r *= v
    return r
def D_Yn(f, x, i, n, th=th):
    f = f*prod(x[:i])**n
    f = D_Y(f, x, i, th)
    return sp.cancel(f*prod(x[:i-1])**(-n))
def D_Yn_inv(f, x, i, n, th=th):  # (Y_{i,n})^-1 = (X_1..X_i)^-n Y_i^-1 (X_1..X_{i-1})^n
    f = f*prod(x[:i-1])**n
    f = D_Yinv(f, x, i, th)
    return sp.cancel(f*prod(x[:i])**(-n))

# ---------- Hikita ----------
def H_T(f, x, i, t):
    s = swap(f, x, i)
    return sp.cancel(t*s + (t-1)*(f-s)/(1-x[i-1]/x[i]))
def H_Ti(f, x, i, t):
    s = swap(f, x, i); r = x[i-1]/x[i]
    return sp.cancel(s + (1-1/t)*(r*f-s)/(1-r))
def H_Pi(f, x):
    return sp.cancel(x[0]*sub(f, x, list(x[1:]) + [x[0]/q]))
def H_Y(f, x, i, t):  # t^{m-i} T_{i-1}..T_1 Pi T_{m-1}^-1 .. T_i^-1
    m = len(x)
    for j in range(i, m): f = H_Ti(f, x, j, t)
    f = H_Pi(f, x)
    for j in range(1, i): f = H_T(f, x, j, t)
    return sp.cancel(t**(m-i)*f)

def monos(x, deg=3, neg=True):
    out = []
    N = len(x)
    for e in itertools.product(range(0, deg+1), repeat=N):
        if sum(e) <= deg: out.append(prod([v**k for v, k in zip(x, e)]))
    if neg:
        out += [x[0]**-1, x[-1]**-2*x[0], x[0]*x[-1]**-1]
    return out

def is_mult(op, x, ms):
    """return (True,g) if op(f)=g f for all f in ms, g=op(1)"""
    g = simp(op(sp.Integer(1)))
    bad = 0; ratios = set()
    for f in ms:
        r = simp(op(f)/f); ratios.add(r)
        if simp(r - g) != 0: bad += 1
    return bad == 0, g, bad, len(ratios)

def checks(N, x):
    t = th**2
    # dictionary checks
    for f in monos(x, 2, False)[:6]:
        for i in range(1, N):
            assert simp(H_T(f, x, i, t) - th*D_T(f, x, i)) == 0
            assert simp(H_Ti(f, x, i, t) - D_Ti(f, x, i)/th) == 0
            assert simp(D_Ti(D_T(f, x, i), x, i) - f) == 0
            assert simp(H_Ti(H_T(f, x, i, t), x, i, t) - f) == 0
        assert simp(H_Pi(f, x) - x[0]*D_pi(f, x)) == 0
        assert simp(D_pi(D_pii(f, x), x) - f) == 0
        for i in range(1, N+1):
            assert simp(D_Yinv(D_Y(f, x, i), x, i) - f) == 0
    # Hikita Y commute (C1), DFK Y_{i,n} commute (C2) on a monomial
    f = x[0]**2*x[-1]
    for i, j in itertools.combinations(range(1, N+1), 2):
        assert simp(H_Y(H_Y(f, x, j, t), x, i, t) - H_Y(H_Y(f, x, i, t), x, j, t)) == 0
        for n in (-1, 0, 1):
            assert simp(D_Yn(D_Yn(f, x, j, n), x, i, n) - D_Yn(D_Yn(f, x, i, n), x, j, n)) == 0
    print(f"N={N}: dictionary T^H=theta*T, (T^H)^-1=theta^-1 T^-1, Pi=x1*pi, inverses, C1,C2 commutation: OK")

def run(N, deg):
    x = mk(N); t = th**2
    checks(N, x)
    ms = monos(x, deg)
    print(f"N={N}: {len(ms)} test monomials (deg<= {deg} plus 3 Laurent)")
    for i in range(1, N+1):
        for n in (-2, -1, 0, 1, 2):
            if i >= 2 and n in (-2, 2) and N > 3: continue
            ok, g, bad, nr = is_mult(lambda f: H_Y(D_Yn(f, x, i, n), x, i, t), x, ms)
            pred = th**(N-1)*x[0]**(n+1) if i == 1 else None
            extra = ''
            if pred is not None: extra = ' matches theta^(N-1) x1^(n+1): %s' % (simp(g - pred) == 0)
            print(f"  Y^H_{i} o Y^DFK_({i},{n}): mult={ok} bad={bad}/{len(ms)} distinct ratios={nr}  op(1)={g if ok else '-'}{extra}")
    # Exact identity Y^H_i = theta^(N-1) (Y_{i,-1})^{-1}, as operators
    for i in range(1, N+1):
        bad = sum(simp(H_Y(f, x, i, t) - th**(N-1)*D_Yn_inv(f, x, i, -1)) != 0 for f in ms)
        print(f"  Y^H_{i} == theta^(N-1) (Y^DFK_({i},-1))^-1 : bad={bad}/{len(ms)}")
        bad = sum(simp(H_Y(f, x, i, t) - th**(N-1)*D_Yn_inv(f, x, i, 1)) != 0 for f in ms)
        print(f"  Y^H_{i} == theta^(N-1) (Y^DFK_({i},+1))^-1 : bad={bad}/{len(ms)}")
    # controls (must break)
    print("  controls (i=1, n=0,1):")
    for n in (0, 1):
        r = is_mult(lambda f: H_Y(D_Yn_inv(f, x, 1, n), x, 1, t), x, ms)
        print(f"   K1 no inversion: Y^H_1 o (Y_(1,{n}))^-1 mult={r[0]} bad={r[2]}")
        r = is_mult(lambda f: H_Y(D_Yn(f, x, 1, n), x, 1, th), x, ms)
        print(f"   K2 t=theta (not theta^2): mult={r[0]} bad={r[2]}")
        r = is_mult(lambda f: H_Y(D_Yn(f, x, 1, n, th=1/th), x, 1, t), x, ms)
        print(f"   K3 DFK theta->theta^-1: mult={r[0]} bad={r[2]}")
        r = is_mult(lambda f: D_Yn(H_Y(f, x, 1, t), x, 1, n), x, ms)
        print(f"   K4 reversed order Y^DFK_(1,{n}) o Y^H_1: mult={r[0]} bad={r[2]} op(1)={r[1] if r[0] else '-'}")

if __name__ == '__main__':
    for N, deg in [(2, 3), (3, 3), (4, 2)]:
        run(N, deg); sys.stdout.flush()
