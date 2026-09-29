"""Test MacBeth's claim (UID 297, sec. 4-5): for fixed kernel brace D, quotient H,
actions (mu, sigma) and residue [nu], EVERY additive class [beta] admits a tau solving (b),(c).

Setting H = Z/2 (trivial brace), D = index-2 ideal, (D,+) = Z/m cyclic, (D,o) abelian.
(G,+) abelian => mu trivial; H^2((Z/2,+), D; mu=1) = D/2D = Z/2:
   beta = 0  <=>  (G,+) = Z/2 x Z/m,   beta != 0  <=>  (G,+) = Z/2m.
An extension with prescribed (D-brace, sigma, [nu]) and beta!=0 exists iff some brace on Z/2m
realises that data.  So: realised-beta-classes per data tuple.
"""
import sys
from math import gcd
from enum_braces import AbGroup, enumerate_braces


def analyse(mods, m):
    G = AbGroup(mods)
    braces, A = enumerate_braces(G)
    n = G.n
    res = []
    # index-2 cyclic subgroups of order m
    def mult(g, c):
        v = G.zero
        for _ in range(c):
            v = G.add[v][g]
        return v
    Dsets = {}
    for g in range(n):
        if G._order(g) == m:
            D = frozenset(mult(g, c) for c in range(m))
            Dsets.setdefault(D, []).append(g)
    for lam in braces:
        L = [A[i] for i in lam]
        def circ(a, b):
            return G.add[a][L[a][b]]
        cinv = {}
        for a in range(n):
            for b in range(n):
                if circ(a, b) == G.zero:
                    cinv[a] = b
        for D, gens in Dsets.items():
            if not all(set(L[a][d] for d in D) == D for a in range(n)):
                continue  # not a left ideal => not an ideal (index 2 => normal in o)
            if any(circ(a, b) != circ(b, a) for a in D for b in D):
                continue  # (D,o) nonabelian
            outside = [k for k in range(n) if k not in D]
            beta_zero = any(G.add[k][k] in {G.add[d][d] for d in D} for k in outside)
            omega_zero = any(G.add[k][k] == G.zero and circ(k, k) == G.zero for k in outside)
            best = None
            for g in gens:  # generator choice = identification D = Z/m
                coord = {mult(g, c): c for c in range(m)}
                def mul_of(aut):  # aut restricted to D as multiplier
                    return coord[aut[g]]
                lamD = tuple(mul_of(L[mult(g, x)]) for x in range(m))
                nu = tuple(sorted(set(mul_of(L[k]) for k in outside)))
                sig = tuple(sorted(set(tuple(coord[circ(cinv[k], circ(mult(g, x), k))]
                                             for x in range(m)) for k in outside)))
                Lsub = tuple(sorted(set(lamD)))
                key = (lamD, Lsub, nu, sig)
                if best is None or key < best:
                    best = key
            res.append((best, beta_zero, omega_zero))
    return res


if __name__ == "__main__":
    m = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    table = {}
    for mods in [(2 * m,), (2, m)]:
        for key, bz, oz in analyse(mods, m):
            table.setdefault(key, set()).add((mods, 'beta=0' if bz else 'beta!=0',
                                              'Omega=0' if oz else 'Omega!=0'))
    print(f"D = Z/{m}; data tuple (lambda^D, L, [nu], sigma) -> realised (G+, beta, Omega)")
    for key in sorted(table):
        lamD, Lsub, nu, sig = key
        betas = {b for _, b, _ in table[key]}
        flag = '' if betas == {'beta=0', 'beta!=0'} else '   <-- beta!=0 NOT realised' if betas == {'beta=0'} else '   <-- beta=0 NOT realised'
        print(f"L={Lsub} [nu]={nu} sigma={sig} lamD={lamD}")
        print("    ", sorted(table[key]), flag)
