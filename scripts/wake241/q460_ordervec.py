"""Order-vectors of printed tie-free Delta along y=x, y=qx, x=qy (Clio Lemma 6) and the q where they change."""
import sympy as sp
x, y, t, q = sp.symbols('x y t q')
A = (y-q**2*x)*(y-q**3*x)/((y-x)*(y-q*x)); B = (y-q**2*x)*(q**2*y-x)/(y-x)**2
C = (y-q**2*x)**2*(x-q*y)/((y-x)**2*(y-q*x)); sw = lambda f: f.subs({x: y, y: x}, simultaneous=True)
F = {'A': A, 'AT': sw(A), 'B': B, 'BT': sw(B), 'C': C, 'CT': sw(C)}
def order(f, sub, qv):
    g = sp.factor(sp.together(f.subs(q, qv).subs(sub)))
    n, d = sp.fraction(g); o = lambda p: min(m[0] for m in sp.Poly(sp.expand(p), t).monoms())
    return o(n) - o(d)
divs = [{y: x+t}, {y: q*x+t}, {x: q*y+t}]
for name, qv in [('generic', q), ('q=-1', -1), ('q^2-q+1=0', sp.Rational(1,2)+sp.sqrt(3)*sp.I/2),
                 ('q prim cube', sp.Rational(-1,2)+sp.sqrt(3)*sp.I/2), ('q=i', sp.I)]:
    vec = {k: tuple(order(f, {kk: vv.subs(q, qv) for kk, vv in d.items()}, qv) for d in divs) for k, f in F.items()}
    print(f'{name:12s}', vec, 'injective on 5 classes:', len(set(vec.values())) == 5)
print('Delta at q=-1:', [sp.simplify(F[k].subs(q, -1)) for k in ('A', 'B', 'C')])
