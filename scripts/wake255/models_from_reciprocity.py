# Wake 255 / Q36: check that the LS-parity model F_kappa(2j) and the OG model Omega(P_kappa, j+1),
# extended as polynomials in j to negative j, agree with the SR_LG / SR_OG strict counts:
#   F_kappa(2j)|_{j=-t}   == (-1)^{|kappa|} #{strict P_kappa -> [2t-1], diag even}     (t >= 1)
#   Omega(P_kappa,j+1)|_{j=-t} == (-1)^{|kappa|} #{strict P_kappa -> [t-1]}             (t >= 1)
# Counts are exact DP over order ideals (= strict rho inside kappa), one value layer at a time.
# Optional: compare also to H_kappa(-t) from Demazure characters (prove249d/dem.py, wake253/og.py).
import sys, itertools
from fractions import Fraction
import sympy as sp
sys.path.insert(0, '../prove249d'); sys.path.insert(0, '../wake253')

def boxes(k): return frozenset((i, i + c) for i, li in enumerate(k) for c in range(li))
def ideals(k):  # order ideals of P_kappa (NW order) = D(rho), rho strict with rho_i <= kappa_i
    out = [()]
    def rec(i, prev, cur):
        for x in range(1, min(k[i], prev - 1) + 1):
            out.append(cur + (x,))
            if i + 1 < len(k): rec(i + 1, x, cur + (x,))
    if k: rec(0, 10**9, ())
    return sorted({boxes(r) for r in out}, key=len)

def count(k, values, diag_parity, strict):
    """# maps f: P_kappa -> values (increasing list), order-preserving (strict or weak),
    with f(diag) of parity diag_parity (None = no condition)."""
    P = boxes(k); J = ideals(k)
    def ok_layer(L, v):
        if diag_parity is not None and v % 2 != diag_parity and any(a == b for a, b in L): return False
        if strict:  # layer must be an antichain
            for (a, b) in L:
                for (c, d) in L:
                    if (a, b) != (c, d) and a <= c and b <= d: return False
        return True
    vec = {frozenset(): 1}
    for v in values:
        nv = {}
        for I, c in vec.items():
            for K in J:
                if I <= K and ok_layer(K - I, v): nv[K] = nv.get(K, 0) + c
        vec = nv
    return vec.get(P, 0)

j = sp.symbols('j')
def interp(pts):
    return sp.expand(sp.interpolate(pts, j))

def main(n, T=6, dem_upto=0):
    bad = 0; tot = 0; polyfail = 0; ctrl = {'LG[2t]': 0, 'LG-noparity': 0, 'LG F(2j+1)': 0, 'OG[t]': 0}
    demtot = 0; dembad = 0
    if dem_upto:
        from dem import chi as chiC
        import og
    for L in range(n + 1):
        for c in itertools.combinations(range(1, n + 1), L):
            k = tuple(sorted(c, reverse=True)); d = sum(k); s = (-1) ** d
            # LG model F_kappa(2j): weak maps into {0..2j}, diag even; fit on j=0..d, test j=d+1,d+2 (polynomiality)
            F = [count(k, range(0, 2 * x + 1), 0, False) for x in range(d + 3)]
            pF = interp([(x, F[x]) for x in range(d + 1)])
            if any(pF.subs(j, x) != F[x] for x in (d + 1, d + 2)): polyfail += 1
            Fodd = [count(k, range(0, 2 * x + 2), 0, False) for x in range(d + 1)]
            pFodd = interp([(x, Fodd[x]) for x in range(d + 1)])
            # OG model Omega(P, j+1): weak maps into {0..j}
            O = [count(k, range(0, x + 1), None, False) for x in range(d + 3)]
            pO = interp([(x, O[x]) for x in range(d + 1)])
            if any(pO.subs(j, x) != O[x] for x in (d + 1, d + 2)): polyfail += 1
            for t in range(1, T + 1):
                NLG = count(k, range(1, 2 * t), 0, True)
                NOG = count(k, range(1, t), None, True)
                tot += 1
                lg = pF.subs(j, -t); ogv = pO.subs(j, -t)
                if lg != s * NLG or ogv != s * NOG:
                    bad += 1; print('FAIL', k, t, lg, s * NLG, ogv, s * NOG)
                if lg != s * count(k, range(1, 2 * t + 1), 0, True): ctrl['LG[2t]'] += 1
                if lg != s * count(k, range(1, 2 * t), None, True): ctrl['LG-noparity'] += 1
                if pFodd.subs(j, -t) != s * NLG: ctrl['LG F(2j+1)'] += 1
                if ogv != s * count(k, range(1, t + 1), None, True): ctrl['OG[t]'] += 1
                if t <= dem_upto:
                    demtot += 1
                    hl = chiC(set(k), n, -t)[0]; ho = og.chi_dem(k, n, -t)
                    if hl != lg or ho != ogv: dembad += 1; print('DEMFAIL', k, t, hl, lg, ho, ogv)
    print(f'parts<={n}: {tot} (kappa,t) pairs, t=1..{T}: model-vs-strict failures {bad}; '
          f'polynomiality failures (2 extra pts beyond fit) {polyfail}')
    print('  negative controls (failures, should be >0):', ctrl)
    if dem_upto: print(f'  Demazure H(-t), t<= {dem_upto}: {demtot} pairs, {dembad} failures')

if __name__ == '__main__':
    n = int(sys.argv[1]); T = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    dem = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    main(n, T, dem)
