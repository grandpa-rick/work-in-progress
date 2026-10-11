# Day 235 printed check of FPSAC Thm 2.2 (Dominance support), re-implemented from the PRINTED text only:
#   Prop 2.1:  E_k F = sum_{|A|=k} c_A X_A F(X_{A^c}, s X_A),  c_A = prod_{i in A, j notin A} (x_i - t x_j)/(x_i - x_j)
#   e*_lam := E_{lam_1} ... E_{lam_l}(1) = sum_mu c_{lam,mu}(s,t) e_mu
#   Thm 2.2 claims: (i) c in Q[s,t] (no t- or s-denominators); (ii) support = {mu ⊵ lam}, c_{lam,lam}=s^{n(lam)};
#   (iii) c|_{s=1}=0 for mu != lam; (iv) c != 0 for every mu ⊵ lam; (v) s-adic valuation of c_{lam,mu} = n(mu);
#   (vi) at t=0, [s^{n(mu)}] c_{lam,mu} = 1.
# Fully symbolic in x_1..x_N, s, t. Coefficient domain is the FIELD Q(s,t), so polynomiality is tested, not assumed.
# No earlier engine is imported for the computation; scripts/day234 table is only used for a final cross-comparison.
import itertools, sys, os, time
from sympy import QQ
from sympy.polys.rings import ring
from sympy.polys.fields import field

K, S, T = field("s,t", QQ)

def partitions(n, mx=None):
    if mx is None: mx = n
    if n == 0: yield (); return
    for k in range(min(n, mx), 0, -1):
        for r in partitions(n - k, k): yield (k,) + r

def conj(l): return tuple(sum(1 for x in l if x > i) for i in range(l[0])) if l else ()
def nfun(l): return sum(i * x for i, x in enumerate(l))
def dom(mu, lam):
    a = b = 0
    for i in range(max(len(mu), len(lam))):
        a += mu[i] if i < len(mu) else 0; b += lam[i] if i < len(lam) else 0
        if a < b: return False
    return sum(mu) == sum(lam)

class Sym:
    def __init__(self, N, variant='printed'):
        self.N = N; self.variant = variant
        self.R, *self.x = ring(",".join(f"x{i}" for i in range(N)), K)
        x = self.x
        self.V = self.R(1)
        for i in range(N):
            for j in range(i + 1, N): self.V *= (x[i] - x[j])
        self.e = [self.R(1)]
        for k in range(1, N + 1):
            self.e.append(sum((self.R(1) * __import__('math').prod(x[i] for i in A) for A in itertools.combinations(range(N), k)), self.R(0)))
    def Vsub(self, idx):
        v = self.R(1)
        for a, b in itertools.combinations(idx, 2): v *= (self.x[a] - self.x[b])
        return v
    def E(self, k, F):
        N, x = self.N, self.x
        tot = self.R(0)
        for A in itertools.combinations(range(N), k):
            Ac = [j for j in range(N) if j not in A]
            sign = (-1) ** sum(1 for i in A for j in Ac if i > j)
            num = self.R(1)
            for i in A:
                for j in Ac: num *= (x[i] - T * x[j])
            XA = self.R(1)
            for i in A: XA *= x[i]
            if self.variant == 'scale_complement':      # NEGATIVE CONTROL 2: s acts on X_{A^c} instead of X_A
                sub = [(x[j], S * x[j]) for j in Ac]
            else:
                sub = [(x[i], S * x[i]) for i in A]
            Fs = F.compose(sub) if sub else F
            tot += sign * num * self.Vsub(A) * self.Vsub(Ac) * XA * Fs
        q, r = tot.div(self.V)
        assert r == 0, "symmetrisation not divisible by Vandermonde"
        if self.variant == 'divide_qint':               # NEGATIVE CONTROL 1: operator rescaled by 1/[k]_t
            q = q * ((1 - T) / (1 - T ** k))
        return q
    def to_e(self, F):
        out = {}
        F = F.copy() if hasattr(F, 'copy') else F
        while F != 0:
            mon = max(F.monoms())          # lex leading monomial
            c = F.coeff(self.R.from_dict({mon: K(1)}))
            assert all(mon[i] >= mon[i + 1] for i in range(len(mon) - 1)), "not symmetric"
            alpha = tuple(m for m in mon if m > 0)
            mu = conj(alpha)
            if mu and mu[0] > self.N: raise ValueError
            g = self.R(1)
            for m in mu: g *= self.e[m]
            out[mu] = out.get(mu, K(0)) + c
            F = F - c * g
        return {m: v for m, v in out.items() if v != 0}
    def estar(self, lam):
        F = self.R(1)
        for k in reversed(lam): F = self.E(k, F)
        return self.to_e(F)

def is_poly(c):  # c in Q(s,t): denominator constant?
    return c.denom.is_ground
def sval(cnum):  # s-adic valuation of polynomial numerator (PolyElement in s,t)
    return min(m[0] for m in cnum.monoms())

def check(lam, C):
    n = sum(lam); fails = []
    allmu = [mu for mu in partitions(n)]
    # (i)
    for mu, c in C.items():
        if not is_poly(c): fails.append(('i:not-in-Q[s,t]', mu, str(c)))
    if fails: return fails        # remaining claims are about polynomials
    P = {mu: c.numer * (1 / c.denom.LC) for mu, c in C.items()}
    # (ii) support and diagonal
    for mu in P:
        if not dom(mu, lam): fails.append(('ii:support', mu))
    diag = P.get(lam)
    if diag is None or diag != (K.ring.gens[0] ** nfun(lam)): fails.append(('ii:diag', str(diag)))
    for mu in allmu:
        c = P.get(mu)
        if dom(mu, lam):
            if c is None: fails.append(('iv:zero', mu)); continue
            if mu != lam and c.evaluate(K.ring.gens[0], 1) != 0: fails.append(('iii:s=1', mu))
            if sval(c) != nfun(mu): fails.append(('v:sval', mu, sval(c), nfun(mu)))
            c0 = c.evaluate(K.ring.gens[1], 0)
            lead = sum((v for m, v in c0.terms() if m[0] == nfun(mu)), QQ(0)) if c0 != 0 else None
            if lead != 1: fails.append(('vi:t=0lead', mu, lead))
    return fails

if __name__ == '__main__':
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    results = {}
    for variant in ('printed', 'divide_qint', 'scale_complement'):
        tot = bad = 0; kinds = {}
        t0 = time.time()
        for n in range(1, nmax + 1):
            Sy = Sym(n, variant)
            for lam in partitions(n):
                C = Sy.estar(lam); tot += 1
                if variant == 'printed': results[lam] = C
                f = check(lam, C)
                if f:
                    bad += 1
                    for x in f: kinds[x[0]] = kinds.get(x[0], 0) + 1
                    if variant == 'printed': print('FAIL', lam, f)
        print(f"[{variant}] |lam|<={nmax}: {tot} partitions, {tot-bad} pass all claims (i)-(vi), {bad} fail; "
              f"failure kinds: {kinds}  ({time.time()-t0:.1f}s)")
    # cross-comparison with the day234 numeric-x engine table (independent evaluation strategy)
    try:
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'day234'))
        import pickle, sympy as sp
        D = pickle.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'day234', 'ctable_n7.pkl'), 'rb'))
        s_, t_ = sp.symbols('s t'); mism = cmpd = 0
        for lam, C in results.items():
            for mu in set(C) | set(D[lam]):
                a = sp.sympify(str(C[mu].as_expr())) if mu in C else 0
                b = sp.sympify(D[lam].get(mu, '0'))
                cmpd += 1; mism += sp.expand(a - b) != 0
        print(f"cross-check vs scripts/day234/ctable_n7.pkl: {cmpd} coefficients, {mism} mismatches")
    except Exception as ex:
        print('cross-check skipped:', ex)
