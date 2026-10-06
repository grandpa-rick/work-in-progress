# Residue/string expansion of R(u,v) = (1/phi_a) CT[g Z S(u)S(v) K],  S(u)=sum_i 1/(u z_i - 1).
# Configurations: words in {U,V}^a. grade: computed
import sympy as sp, itertools
from sym2pt import Phi, phia
t, u, v = sp.symbols('t u v')
def gval(nu, pts):
    return sp.prod([sum(p**r for p in pts) for r in nu]) if nu else sp.Integer(1)
def config_weight(word, nu, onestring_src=None):
    a = len(word); pts = [None]*a; cnt = {'U':0, 'V':0}; base = {'U':1/u, 'V':1/v}
    for m, L in enumerate(word):
        pts[m] = base[L]*t**cnt[L]; cnt[L] += 1
    W = sp.Integer(1)
    last = {}
    for m, L in enumerate(word):
        for i in range(m):
            if word[i] == L and pts[m] == t*pts[i]:
                W *= (pts[m]-pts[i])   # successor: residue of 1/(z_m - t z_i) part
            else:
                W *= (pts[m]-pts[i])/(pts[m]-t*pts[i])
    # seeds: residue of 1/(u z - 1) at z=1/u is 1/u
    if 'U' in word: W *= 1/u
    if 'V' in word: W *= 1/v
    if onestring_src is not None:  # extra non-residue S factor at position j
        L, j = onestring_src
        W *= 1/((v if L=='U' else u)*pts[j]-1)
    return W*gval(nu, pts)
def R(a, nu):
    tot = 0
    for word in itertools.product('UV', repeat=a):
        if 'U' in word and 'V' in word:
            tot += config_weight(word, nu)
    for j in range(a):
        tot += config_weight('U'*a, nu, ('U', j)) + config_weight('V'*a, nu, ('V', j))
    return sp.cancel(tot/phia(a))
if __name__ == '__main__':
    for a, nu in [(2,()),(2,(1,)),(2,(2,)),(3,(1,)),(3,(2,1)),(2,(1,1)),(4,(2,))]:
        n = a + sum(nu)
        r = sp.expand(sp.cancel(R(a, nu)*u**n*v**n))  # should be polynomial: sum Phi/(..) u^{n-x} v^{n-y}
        ok = True
        for x in range(1, n):
            y = n-x
            c = sp.Poly(r, u, v).coeff_monomial(u**(n-x)*v**(n-y))
            if sp.simplify(c*(1-t**x)*(1-t**y) - Phi(a, nu, x, y)) != 0: ok = False; print('mismatch', a, nu, x, y)
        print(a, nu, 'string expansion matches Phi:', ok, '| polynomial:', sp.Poly(r, u, v).is_polynomial if hasattr(sp.Poly(r,u,v),'is_polynomial') else True)
