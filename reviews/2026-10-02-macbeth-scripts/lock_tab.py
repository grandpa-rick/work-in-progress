"""Day 217: independent tabulation of MacBeth's three-way lock
   (lambda^D_x - 1) b = 2 delta_1(x),  4b = 2(nu sigmahat - 1) mod 8
for H = Z/2 extensions 0 -> D -> G -> Z/2 -> 0, (G,+) abelian, D cyclic ideal.
Conventions: u = s(1) (section, any u not in D), b = 2u in D,
nu = lambda_u|_D, sigmahat(x) = u^{-1} o x o u, delta(x) = lambda_x(u) - u.
"""
import sys
from collections import defaultdict
from braces import AbGroup, regular_subgroups, brace_ops

def cyclic_subgroups(A, order):
    subs = {}
    for g in range(A.n):
        if A.order(g) == order:
            S = frozenset(A.mul(k, g) for k in range(order))
            subs.setdefault(S, g)
    return subs  # set -> a generator

def coord(A, g, order):
    return {A.mul(k, g): k for k in range(order)}

def as_unit(f_on_D, A, g, order, co):
    """if f restricted to D is multiplication by a unit m (wrt generator g), return m else None"""
    m = co[f_on_D(g)]
    for k in range(order):
        x = A.mul(k, g)
        if co[f_on_D(x)] != (m * k) % order:
            return None
    return m

def run(mods, dorder, verbose=False):
    A = AbGroup(mods)
    R = regular_subgroups(A)
    rows = []
    for bi, lam in enumerate(R):
        circ, inv = brace_ops(A, lam)
        for D, g in cyclic_subgroups(A, dorder).items():
            if A.n // dorder != 2:
                continue
            if not all(lam[a][d] in D for a in range(A.n) for d in D):
                continue  # not lambda-invariant => not an ideal (index 2 => normal in o)
            co = coord(A, g, dorder)
            L = tuple(sorted({as_unit(lambda x, d=d: lam[d][x], A, g, dorder, co) for d in D}))
            for u in range(A.n):
                if u in D:
                    continue
                b = A.add[u][u]; bk = co[b]
                nu = as_unit(lambda x: lam[u][x], A, g, dorder, co)
                sh = as_unit(lambda x: circ[circ[inv[u]][x]][u], A, g, dorder, co)
                delta = lambda x: A.add[lam[x][u]][A.neg[u]]
                nsh = as_unit(lambda x: A.add[x][delta(x)], A, g, dorder, co)
                # lock for every x in D: (lam_x - 1) b == 2 delta(x)
                lock_all = all(A.add[lam[x][b]][A.neg[b]] == A.add[delta(x)][delta(x)] for x in D)
                c0 = co[A.add[g][delta(g)]]  # nu sigmahat evaluated at generator
                lam_g = co[lam[g][g]]
                lock_gen = ((lam_g - 1) * bk - 2 * (c0 - 1)) % dorder == 0
                rows.append(dict(brace=bi, D=D, u=u, L=L, b=bk, nu=nu, sh=sh, nsh=nsh, c0=c0,
                                 lock_all=lock_all, lock_gen=lock_gen, lam_g=lam_g,
                                 nsh_eq=(nu is not None and sh is not None and nsh is not None
                                         and (nu * sh) % dorder == nsh)))
    return A, R, rows

if __name__ == "__main__":
    for mods, dorder in [((16,), 8), ((2, 8), 8), ((8,), 4), ((2, 4), 4)]:
        A, R, rows = run(mods, dorder)
        print(f"=== (G,+)=Z{mods}, |D|={dorder}: {len(R)} braces (labelled regular subgroups)")
        byL = defaultdict(list)
        for r in rows:
            byL[r['L']].append(r)
        for L, rs in sorted(byL.items()):
            npairs = len({(r['brace'], r['D']) for r in rs})
            nbr = len({r['brace'] for r in rs})
            print(f"  L={L}: braces={nbr} (brace,D) pairs={npairs} (brace,D,u) triples={len(rs)}")
            print(f"     lock_all {sum(r['lock_all'] for r in rs)}/{len(rs)}, lock_at_gen {sum(r['lock_gen'] for r in rs)}/{len(rs)},"
                  f" nu unit {sum(r['nu'] is not None for r in rs)}, sigmahat unit {sum(r['sh'] is not None for r in rs)},"
                  f" nu*sh additive {sum(r['nsh'] is not None for r in rs)}, nsh=nu*sh {sum(r['nsh_eq'] for r in rs)}")
            # parity of b vs coset of nsh in Aut/L
            tab = defaultdict(lambda: defaultdict(int))
            for r in rs:
                coset = tuple(sorted({(r['nsh'] * l) % dorder for l in L})) if r['nsh'] else None
                tab[coset][r['b'] % 2] += 1
            for coset, d in sorted(tab.items(), key=str):
                print(f"     nsh-coset {coset}: b even {d[0]}, b odd {d[1]}, offset 2(c0-1) values "
                      f"{sorted({(2*(r['c0']-1)) % dorder for r in rs if r['nsh'] and tuple(sorted({(r['nsh']*l)%dorder for l in L}))==coset})}")
