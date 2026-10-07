"""Checks for review of MacBeth 5f1531d (container-derivative uniqueness + tangent notes)."""
import itertools, sympy as sp

# (1) Base change vs dual-number substitution on a polynomial p in N[y].
y, v, e = sp.symbols('y v epsilon')
def trunc(expr):  # mod eps^2
    expr = sp.expand(expr); return sp.expand(expr.subs(e**2, 0).subs(e**3, 0).subs(e**4,0))
for p in [y**2, y**3, y**2 + 1]:
    sub = trunc(p.subs(y, y + e*v))
    print("p =", p, "| dual-number subst:", sub, "| base change N[y]->W[y] image:", p,
          "| p (x) 1 in N[y](x)W viewed as module elt: p*1 (eps-part 0)")

# (2) Ranks in finite free N-modules: T=(-)(x)D, D=N+M, rank M=n, M.M=0.
# D(x)D basis: 1, e_i (inner), f_j (outer), e_i f_j.  T(p)=p(x)D applied to INNER factor kills e_i, e_i f_j.
# V = preimage of zero-section image (span{1}) = span{1, e_i, e_i f_j}? -> depends on which factor;
# compute both and the fibre product D_2 = D x_N D.
def ranks(n):
    basis = ['1'] + [f'e{i}' for i in range(n)] + [f'f{j}' for j in range(n)] + [f'e{i}f{j}' for i in range(n) for j in range(n)]
    def pT(b):  # p (x) D : kill inner e
        if b == '1': return '1'
        if b.startswith('e'): return None   # -> 0
        return 'eps'                          # f_j -> e_j in D
    # over N, x lands in span{1} iff its coefficients on basis elts mapping to non-1 nonzero vanish
    V = [b for b in basis if pT(b) in ('1', None)]
    D2 = 1 + 2*n
    return len(basis), D2, len(V)
for n in range(6):
    tot, D2, V = ranks(n)
    print(f"n={n}: rank D(x)D={tot}=(1+n)^2? {tot==(1+n)**2}; rank D2={D2}; rank V={V}; iso? {D2==V}")

# (3) Poly-morphism model (linear containers, morphisms = functions on shapes): p-hat sends all shapes to the one shape.
for n in range(6):
    shapes = list(itertools.product(range(1+n), repeat=2))
    Pp = [s for s in shapes if s[1] == 0]  # preimage of zero shape under (a,b)->b ... but p-hat o pi: every shape maps to the single y-shape
    print(f"n={n}: Poly-function model |D x_y D|={(1+n)**2}, |P'|={1+n}")

# (4) Flip axiom c.l = l with l: M -> M(x)M an iso of free N-modules.
# Isos of free N-modules permute bases (the atoms). l(e_k)=e_i(x)e_j; c.l=l forces i=j; surjectivity forces rank<=1.
def flip_ok(n):
    targets = [(i,j) for i in range(n) for j in range(n)]
    if len(targets) != n: return False  # no basis bijection at all
    for perm in itertools.permutations(targets):
        if all(i == j for (i, j) in perm): return True
    return False
print("flip-compatible iso l:M->M(x)M exists for rank n:", {n: flip_ok(n) for n in range(5)})

# (5) Canonical v : T_2 -> V for D=W over N-modules, A = N^2, on a box: v(a,b,c) = a + eps_in c + eps1eps2 l(b).
# v is a bijection onto V iff l is; with l = id (rank 1) check on a box.
box = range(4)
imgs = set()
for a in itertools.product(box, repeat=2):
    for b in itertools.product(box, repeat=2):
        for c in itertools.product(box, repeat=2):
            imgs.add((a, c, b))   # components (1, eps_in, eps1eps2); eps_out-component 0
print("D=W, A=N^2: v injective on box:", len(imgs) == 4**6)

# (6) Remark 6 arithmetic: (y+1)<|(y+1)=y+2 ; Dirichlet (y+1)(x)(y+1)=y+3.
print("(y+1)<|(y+1) =", sp.expand((y+1).subs(y, y+1)))
# Dirichlet: y^A (x) y^B = y^{AxB}; exponents 1,0 -> products 1,0,0,0
print("Dirichlet:", sum(y**(a*b) for a in (1,0) for b in (1,0)))

# (7) CD.6 / CD.7 forms as stated in tangent note, checked on generic polynomials.
x, xp, vv, vp, a, b, c, d = sp.symbols('x xp vv vp a b c d')
for f in [x**2, x**3, x**4 + 2*x + 3]:
    Df = lambda X, V: sp.diff(f, x).subs(x, X) * V
    DDf = lambda X, V, Xp, Vp: sp.expand(sp.diff(Df(x, vv), x).subs({x: X, vv: V}) * Xp + sp.diff(Df(x, vv), vv).subs({x: X, vv: V}) * Vp)
    cd6 = sp.simplify(DDf(a, b, 0, d) - Df(a, d)) == 0
    cd7 = sp.simplify(DDf(a, b, c, 0) - DDf(a, c, b, 0)) == 0
    print("f =", f, "CD.6 (as stated, b free):", cd6, " CD.7:", cd7)

# (8) Prop 9: linear iff N(f)=a1 x : solve sum a_k k v x^{k-1} = sum a_k v^k
A = sp.symbols('a0:5')
f = sum(A[k]*x**k for k in range(5))
eq = sp.Poly(sp.expand(sp.diff(f, x)*vv - f.subs(x, vv)), x, vv)
print("Prop 9 constraints:", sp.solve(eq.coeffs(), A, dict=True))
