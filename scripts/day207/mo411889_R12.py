"""Day 207: MO 411889 check. Straighten R_(1,2) (Stanley's comment:
R_(1,2) = (t^3-1) R_(1,1,1) + t R_(2,1)) and Jing's Q_(1,2).

R_a(x;t) = sum_{w in S_n} w( x^a prod_{i<j} (x_i - t x_j)/(x_i - x_j) )   (Macdonald III.1)
P_l = R_l / v_l(t),  v_l = prod_{i>=0} v_{m_i},  v_m = prod_{j=1}^m (1-t^j)/(1-t)  (m_0 = n - l(l))
Q_l = b_l(t) P_l,  b_l = prod_{i>=1} prod_{j=1}^{m_i} (1-t^j)
Jing: Q_a = [z^a] prod_i Q(z_i) prod_{i<j} (1-z_j/z_i)/(1-t z_j/z_i),  Q(z)=prod_k (1-t x_k z)/(1-x_k z)
"""
import itertools, sympy as sp
t = sp.symbols('t')

def R(a, xs):
    n = len(xs); a = tuple(a) + (0,)*(n-len(a))
    tot = 0
    for p in itertools.permutations(range(n)):
        y = [xs[i] for i in p]
        term = sp.Mul(*[y[i]**a[i] for i in range(n)])
        for i in range(n):
            for j in range(i+1, n):
                term *= (y[i]-t*y[j])/(y[i]-y[j])
        tot += term
    return sp.expand(sp.cancel(sp.together(tot)))

def vm(m): return sp.prod([(1-t**j)/(1-t) for j in range(1, m+1)])
def mult(l, n):
    l = [x for x in l if x > 0]
    return {k: l.count(k) for k in set(l)}, n-len(l)
def P(l, xs):
    ms, m0 = mult(l, len(xs))
    return sp.expand(sp.cancel(R(l, xs)/(sp.prod([vm(m) for m in ms.values()])*vm(m0))))
def b(l):
    ms, _ = mult(l, 0)
    return sp.prod([sp.prod([1-t**j for j in range(1, m+1)]) for m in ms.values()])
def Q(l, xs): return sp.expand(b(l)*P(l, xs))

def qn(N, xs):  # q_n = coefficient of z^n in prod (1-t x z)/(1-x z)
    z = sp.symbols('z')
    g = sp.prod([(1-t*x*z)*sum((x*z)**k for k in range(N+1)) for x in xs])
    return [sp.expand(sp.expand(g).coeff(z, k)) for k in range(N+1)]

def jing2(a1, a2, xs):  # [z^a1 w^a2] Q(z)Q(w)(1-w/z)/(1-tw/z) = sum_k f_k q_{a1+k} q_{a2-k}
    q = qn(a1+a2, xs)
    f = lambda k: 1 if k == 0 else t**k - t**(k-1)
    return sp.expand(sum(f(k)*q[a1+k]*q[a2-k] for k in range(a2+1)))

def straighten(F, basis, xs):
    cs = sp.symbols('c0:%d' % len(basis))
    eq = sp.Poly(sp.expand(F - sum(c*B for c, B in zip(cs, basis))), *xs).coeffs()
    sol = sp.solve(eq, cs, dict=True)[0]
    return [sp.factor(sol[c]) for c in cs]

parts3 = [(3,), (2,1), (1,1,1)]
for n in (2, 3, 4):
    xs = sp.symbols('x1:%d' % (n+1))
    print(f"=== n = {n} variables ===")
    for a in [(1,2)] + ([(1,2,0)] if n >= 3 else []) + ([(0,1,2)] if n >= 3 else []):
        if len(a) > n: continue
        Ra = R(a, xs)
        ps = [l for l in parts3 if len(l) <= n]
        c = straighten(Ra, [R(l, xs) for l in ps], xs)
        print(f"R_{a} =", " + ".join(f"({ci})*R_{l}" for ci, l in zip(c, ps)))
    ps = [l for l in parts3 if len(l) <= n]
    J12 = jing2(1, 2, xs)
    print("Jing Q_(2,1) == Q_(2,1):", sp.expand(jing2(2, 1, xs) - Q((2,1), xs)) == 0)
    c = straighten(J12, [Q(l, xs) for l in ps], xs)
    print("Jing Q_(1,2) =", " + ".join(f"({ci})*Q_{l}" for ci, l in zip(c, ps)))

# --- all rearrangements of (1,2,0) in n=3, to locate Stanley's (t^3-1)R_(111)+tR_(21) ---
if __name__ == "__main__":
    xs = sp.symbols('x1:4'); ps = parts3; basis = [R(l, xs) for l in ps]
    print("=== n=3, all rearrangements of (2,1,0) ===")
    for a in sorted(set(itertools.permutations((2,1,0)))):
        c = straighten(R(a, xs), basis, xs)
        print(f"R_{a} =", " + ".join(f"({ci})*R_{l}" for ci, l in zip(c, ps)))
