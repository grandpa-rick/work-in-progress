"""Day 224 closed formulas (symbolic s,t), all from PROVED inputs: Thm 1.5 (Day 223), Column Lemma, recursion (S).
Lam(k,J) = lin_e E_k(e_J) = (-1)^d [k+d]/[k] prod_i pi_k(j_i),  pi_k(j) = [u^j] (-su;t)_k/(-u;t)_k.
star2(k,r) = e_k * e_r (Hikita) = sum_y s^y Lam(k-y,(r-y)) e_{k+r-y} e_y   (Column Lemma + Thm 1.5).
cN1(lam) = c_{lam,(n-1,1)} for l(lam)=3 via (S)."""
import sympy
from functools import lru_cache
s, t = sympy.symbols('s t')
def qi(m): return sympy.Integer(0) + sum(t**i for i in range(m))
@lru_cache(None)
def qbin(m, i):
    if i < 0 or i > m: return sympy.Integer(0)
    num = sympy.Integer(1); den = sympy.Integer(1)
    for j in range(i):
        num *= (1 - t**(m-j)); den *= (1 - t**(j+1))
    return sympy.cancel(num/den)
@lru_cache(None)
def pi(k, j):
    # (-su;t)_k = sum_i s^i t^{C(i,2)} [k,i] u^i ; 1/(-u;t)_k = sum_m (-u)^m [k+m-1, m]
    return sympy.expand(sum(s**i * t**(i*(i-1)//2) * qbin(k, i) * (-1)**(j-i) * qbin(k+j-i-1, j-i) for i in range(0, min(j, k)+1)))
def Lam(k, J):
    J = [j for j in J if j > 0]
    if k == 0: return sympy.Integer(1) if len(J) == 1 else sympy.Integer(0)
    d = sum(J); v = (-1)**d * qi(k+d) / qi(k)
    for j in J: v *= pi(k, j)
    return v
def key(*parts): return tuple(sorted([p for p in parts if p > 0], reverse=True))
@lru_cache(None)
def star2(k, r):
    out = {}
    for y in range(0, min(k, r)+1):
        a, b = k-y, r-y
        C = Lam(a, [b]) if (a > 0 and b > 0) else sympy.Integer(1)
        kk = key(k+r-y, y); out[kk] = out.get(kk, 0) + s**y * C
    return out
def xtop(F, m):
    """[x_1^m] of a dict in e-basis (keys partitions) -> dict in x' e-basis. e_k = e_k' + x_1 e_{k-1}'."""
    out = {}
    import itertools
    for kk, c in F.items():
        L = len(kk)
        for S in itertools.combinations(range(L), m):
            parts = [kk[i]-1 if i in S else kk[i] for i in range(L)]
            if min(parts, default=0) < 0: continue
            k2 = key(*parts); out[k2] = out.get(k2, 0) + c
    return out
def linE(k, G):  # lin_e E_k(G), G dict
    return sum(c * Lam(k, list(J)) for J, c in G.items())
def cN1(lam):
    l1, l2, l3 = lam; n = sum(lam); n1 = l1-1
    F = star2(l2, l3); F1 = xtop(F, 1)
    F2 = {kk: s*c for kk, c in star2(l2-1, l3-1).items()}
    assert all(sympy.simplify(xtop(F, 2).get(kk, 0) - c) == 0 for kk, c in F2.items()), 'column lemma check'
    R1 = s * linE(n1, F1)
    R2 = t**l1 * linE(l1, F2)
    if n1 >= 1:
        R3 = s**2*(1-t) * sum(c * (-1)**sum(J) * qi(n1+sum(J)+1) * sympy.Mul(*[pi(n1, j) for j in J if j > 0]) for J, c in F2.items())
    else:  # K_0(G) = (s e_1 G - e_1 G)/(s-1) = e_1 G ; lin_e(e_1 G) = [G const]
        R3 = s**2*(1-t) * F2.get((), 0)
    cnu = linE(n1, star2(l2-1, l3-1)) if n-3 > 0 else 1
    mult = 3 if n == 3 else 1   # rho(e_{(0,0,0)}) = 3 e_1
    return sympy.cancel(R1 + R2 + R3 - mult * s**3 * cnu)
