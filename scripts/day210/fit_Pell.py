"""Day 210: fit / test prediction P-ell on the symbolic recursion output.
For each term (lam=(b,), I) with coefficient R: divide by c^b prod_{c<c'} K_{i_c i_c'}(z_c,z_c'), c = s^ell t^{-sum I}.
P-ell (shape) says the quotient is a Laurent polynomial sum_n  mono(s,t) * prod_c N^{(n_c)}_{i_c} z_c^{-n_c}, sum n = k-b.
We report, for each monomial, the ratio to prod N (should be a monomial s^alpha t^beta) and compare with
  version A (connection note): s^{sum(n_c-i_c)} t^{-sum_{c!=c'}(n_c-i_c) i_c'}
  version B (GF-natural, X_c = c x/(s z_c)):  s^{(ell-1) sum n_c - sum i_c} t^{-sum_{c!=c'}(n_c-i_c) i_c'}   [= A at ell=2]
"""
import sympy as sp, pickle, sys, itertools
s, t = sp.symbols('s t')
def poch(a, N): return sp.Mul(*[1 - a*t**i for i in range(N)])
def alpha(J): return sp.Integer(0) if J < 0 else sp.Mul(*[s - t**i for i in range(1, J+1)])/poch(t, J)
def cc(n, j):
    if j < 0 or j > n: return sp.Integer(0)
    return poch(s, n-j)/poch(t, n-j)*(alpha(j) - s*t**(n-j)*alpha(j-1))
def Nn(n, j): return t**(-n*j)*cc(n, j)
def K(i, j, z, w):
    return sp.Mul(*[(t**p*z - s*w)/(t**p*z - t**j*w) for p in range(i)]) * sp.Mul(*[(s*z - t**r*w)/(t**i*z - t**r*w) for r in range(j)])

def analyse(ell, k, G, verbose=True):
    Z = sp.symbols(f'z1:{ell+1}')
    allok = {'shape': True, 'A': True, 'B': True}
    nonchain = [key for key in G if len(key[0]) > 1]
    for (lam, I), R in sorted(G.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        if len(lam) > 1: continue
        b = lam[0] if lam else 0
        cst = s**ell*t**(-sum(I))
        Kp = sp.Mul(*[K(I[a], I[bb], Z[a], Z[bb]) for a in range(ell) for bb in range(a+1, ell)])
        Q = sp.factor(sp.cancel(R/(cst**b*Kp)))
        num, den = sp.fraction(sp.cancel(Q))
        # den must be monomial in Z times something in s,t only
        dpoly = sp.Poly(den, *Z)
        if len(dpoly.terms()) != 1:
            allok['shape'] = False
            print(f'  NON-LAURENT  I={I} b={b}: quotient = {Q}'); continue
        (dmon, dco), = dpoly.terms()
        npoly = sp.Poly(sp.expand(num), *Z)
        line = []
        for mon, co in npoly.terms():
            n = [dmon[c] - mon[c] for c in range(ell)]
            if any(x < 0 for x in n) or sum(n) != k - b:
                allok['shape'] = False; line.append(f'BADEXP n={n}'); continue
            pN = sp.Mul(*[Nn(n[c], I[c]) for c in range(ell)])
            if pN == 0:
                allok['shape'] = False; line.append(f'n={n}: prodN=0 but coeff {sp.factor(co/dco)}'); continue
            r = sp.factor(sp.cancel(co/dco/pN))
            tA = -sum((n[c]-I[c])*I[d] for c in range(ell) for d in range(ell) if d != c)
            predA = s**sum(n[c]-I[c] for c in range(ell))*t**tA
            predB = s**((ell-1)*sum(n) - sum(I))*t**tA
            okA = sp.cancel(r - predA) == 0; okB = sp.cancel(r - predB) == 0
            allok['A'] &= okA; allok['B'] &= okB
            if not (okA or okB):
                allok['shape'] &= (sp.Poly(sp.numer(r), s, t).length() == 1)  # still a monomial?
            line.append(f'n={n}: {r} [A:{okA} B:{okB}]')
        # also check all needed monomials appear (missing ones would be zero coefficient where prod N != 0)
        if verbose: print(f'  I={I} b={b}: ' + '; '.join(line))
    # completeness: predicted terms with no recursion term
    for b in range(k+1):
        for I in itertools.product(range(k-b+1), repeat=ell):
            if sum(I) > k-b: continue
            if ((b,) if b else (), I) not in G:
                # predicted coefficient nonzero?
                nz = any(sp.Mul(*[Nn(n[c], I[c]) for c in range(ell)]) != 0
                         for n in itertools.product(range(k-b+1), repeat=ell) if sum(n) == k-b)
                if nz: allok['shape'] = False; print(f'  MISSING predicted term I={I} b={b}')
    print(f'ell={ell} k={k}: non-chain terms: {nonchain if nonchain else "none"};  pairwise-shape: {allok["shape"]};  version A exact: {allok["A"]};  version B exact: {allok["B"]}')
    return allok

if __name__ == '__main__':
    path, ell = sys.argv[1], int(sys.argv[2])
    res = pickle.load(open(path, 'rb'))
    for k in sorted(res):
        if k == 0: continue
        analyse(ell, k, res[k])

# ---- version C (fitted at ell=3,k<=2; = TC at ell=2): s^{(ell-1) sum(n_c-i_c)} t^{-sum_{c!=c'}(n_c-i_c) i_c'} ----
def predC(ell, k, Z):
    """full predicted dict (lam, I) -> coefficient, P-ell with pairwise K and version-C prefactor, c^b dressing."""
    out = {}
    for b in range(k+1):
        for n in itertools.product(range(k-b+1), repeat=ell):
            if sum(n) != k-b: continue
            for I in itertools.product(*[range(nc+1) for nc in n]):
                cst = s**ell*t**(-sum(I))
                tA = -sum((n[c]-I[c])*I[d] for c in range(ell) for d in range(ell) if d != c)
                wt = cst**b * s**((ell-1)*sum(n[c]-I[c] for c in range(ell)))*t**tA
                Kp = sp.Mul(*[K(I[a], I[bb], Z[a], Z[bb]) for a in range(ell) for bb in range(a+1, ell)])
                term = wt*Kp*sp.Mul(*[Nn(n[c], I[c])*Z[c]**(-n[c]) for c in range(ell)])
                key = ((b,) if b else (), tuple(I))
                out[key] = out.get(key, 0) + term
    return out

def exact_compare(ell, k, G):
    Z = sp.symbols(f'z1:{ell+1}')
    P = predC(ell, k, Z)
    bad = [key for key in set(P) | set(G) if sp.cancel(sp.together(P.get(key, 0) - G.get(key, 0))) != 0]
    print(f'ell={ell} k={k}: exact term-by-term compare with P-ell(version C): {"ALL EQUAL" if not bad else "MISMATCH " + str(bad)}  ({len(G)} terms)')
    return not bad
