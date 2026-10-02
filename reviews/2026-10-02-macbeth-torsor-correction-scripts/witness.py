"""Extract explicit witnesses for the crit-but-unbased triplet classes at |G|=16, D=Z/2xZ/4."""
from collections import defaultdict
from braces import AbGroup, regular_subgroups, brace_ops
from criterion import run
groups = [(16,), (2, 8), (4, 4), (2, 2, 4)]
cls = defaultdict(list); data = {}
for mods in groups:
    A, R, out = run(mods); data[mods] = (A, R)
    for r in out: cls[r['inv']].append((mods, r))
for key, mem in cls.items():
    if key[0] != (2, 4): continue
    if any(r['based'] for _, r in mem) or not any(r['crit'] for _, r in mem): continue
    mods, r = mem[0]
    A, R = data[mods]; lam = R[r['brace']]; circ, inv = brace_ops(A, lam)
    D = sorted(r['D']); E = A.elts
    print(f"class |L|={r['Lsize']}, members {len(mem)} all in {sorted({m for m,_ in mem})}")
    print("  D =", [E[x] for x in D])
    gens = [A.idx[(1,0)], A.idx[(0,1)]]
    for gi in gens:
        print(f"  lambda_{E[gi]}: ", {E[x]: E[lam[gi][x]] for x in [A.idx[(1,0)], A.idx[(0,1)]]})
    twoD = {A.add[x][x] for x in D}
    for u in range(A.n):
        if u in r['D']: continue
        b = A.add[u][u]
        good = all(A.add[A.add[lam[d][u]][A.neg[u]]][A.add[lam[d][u]][A.neg[u]]] == A.zero for d in D)
        if good:
            print(f"  good section u={E[u]}: b=2u={E[b]}, 2D={[E[x] for x in twoD]}, delta_1(d) for d in D:",
                  {E[d]: E[A.add[lam[d][u]][A.neg[u]]] for d in D})
            break
    print("  brace check: lam hom", all(lam[circ[a][c]] == tuple(lam[a][lam[c][x]] for x in range(A.n)) for a in range(A.n) for c in range(A.n)))
