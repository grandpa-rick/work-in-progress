"""Print explicit Z/16 braces with D=2Z/16 ~ (Z/8, lambda_x=5^x), sigma=id, [nu]={3,7}; and Z2xZ8 ones with sigma=id."""
from enum_braces import AbGroup, enumerate_braces
for mods in [(16,), (2, 8)]:
    G = AbGroup(mods); br, A = enumerate_braces(G)
    for lam in br:
        L = [A[i] for i in lam]
        circ = lambda a, b: G.add[a][L[a][b]]
        for g in range(G.n):
            if G._order(g) != 8: continue
            Dl = [G.zero]
            for _ in range(7): Dl.append(G.add[Dl[-1]][g])
            D = set(Dl); coord = {d: i for i, d in enumerate(Dl)}
            if not all(set(L[a][d] for d in D) == D for a in range(G.n)): continue
            lamD = tuple(coord[L[d][g]] for d in Dl)
            if lamD != (1, 5, 1, 5, 1, 5, 1, 5): continue
            out = [k for k in range(G.n) if k not in D]
            inv = {a: next(b for b in range(G.n) if circ(a, b) == G.zero) for a in range(G.n)}
            sig = {tuple(coord[circ(inv[k], circ(d, k))] for d in Dl) for k in out}
            nu = sorted({coord[L[k][g]] for k in out})
            if sig == {tuple(range(8))}:
                # lambda as table of images of generators
                gens = [G.idx[tuple(1 if j == i else 0 for j in range(len(mods)))] for i in range(len(mods))]
                tab = {G.elts[a]: tuple(G.elts[L[a][x]] for x in gens) for a in range(G.n)}
                print(mods, "D gen", G.elts[g], "nu", nu)
                print("   lambda_a(generators):", tab)
                break
