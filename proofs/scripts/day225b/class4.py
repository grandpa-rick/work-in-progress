# Day 225: class-4 leads from CLOSED formulas only: 2-point formula (formula2pt), Thm 1.5, Thm 7.1 (Mthm) + derivation.
# [(s-1)^2] c_{lam,(x,y)} = [e_x e_y](Gamma_a(e_b,e_c) + D_a D_b e_c)  (Prop 5.1). grade: computed
import sys; sys.path.insert(0, '../day224')
import sympy as sp
from formula2pt import Phi_formula, t
def qi(m): return sum(t**i for i in range(m))
def pspec(nu, a): return sp.prod([sum(t**(r*i) for i in range(a)) for r in nu]) if nu else sp.Integer(1)
def lin_T(a, nu):  # Thm 1.5: lin_e T_a(p_nu)
    d = sum(nu); n = a + d
    return (-1)**d*qi(n)/qi(a)*pspec(nu, a)
def U(a, nu, x, y):  # [e_x e_y] T_a(p_nu), x+y = a+|nu|
    n = x + y; m = 2 if x == y else 1
    return sp.cancel(((-1)**n*Phi_formula(a, tuple(sorted(nu, reverse=True)), x, y) - lin_T(a, nu))/m)
def L(a, b): return (1-t**(a+b))*(t**(a*b)-1)/((1-t**a)*(1-t**b))
def M(k, r):  # Thm 7.1: D_k(e_r) = B(e_k,e_r) in e-basis dict
    if r == 0 or k == 0: return {}
    if k < r: k, r = r, k
    out = {(k, r): sp.Integer(r)}
    for j in range(1, r+1):
        key = tuple(z for z in sorted((k+j, r-j), reverse=True) if z)
        out[key] = out.get(key, 0) + L(k-r+j, j)
    return out
def Dder(a, F):
    out = {}
    for key, v in F.items():
        for i, z in enumerate(key):
            rest = key[:i] + key[i+1:]
            for k2, w in M(a, z).items():
                kk = tuple(sorted(rest + k2, reverse=True)); out[kk] = out.get(kk, 0) + v*w
    return out
def lead(lam, mu, a_idx=0):
    lam = list(lam); a = lam.pop(a_idx); b, c = lam
    x, y = mu; key = tuple(sorted((x, y), reverse=True))
    tot = (-1)**(b+c)*U(a, (b, c), x, y)       # r=b, q=c
    for r in range(1, b):                        # r<b, q=c : e_{b-r} * lin_e T_a(p_r p_c)
        if tuple(sorted((b-r, a+r+c), reverse=True)) == key:
            tot += (-1)**(r+c)*lin_T(a, (r, c))
    for q in range(1, c):
        if tuple(sorted((c-q, a+b+q), reverse=True)) == key:
            tot += (-1)**(b+q)*lin_T(a, (b, q))
    tot += Dder(a, M(b, c)).get(key, 0)
    return sp.factor(sp.cancel(tot))
if __name__ == '__main__':
    target = 2*t**13+3*t**12+3*t**11+6*t**10+6*t**9+6*t**8+9*t**7+6*t**6+3*t**5+9*t**4+6*t**3+3*t+4
    v = lead((3,3,3), (7,2))
    print('KILL TEST (3,3,3)->(7,2):', sp.expand(v), '  MATCH:', sp.expand(v - target) == 0)
