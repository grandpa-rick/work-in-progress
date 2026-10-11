# Day 234 printed-statement checks, part A (statements about the coefficients c_{lam,mu}(s,t)).
# Ground truth: ctable_n7.pkl from engine.py (subset formula, N=7). Formulas are transcribed from printed.txt.
from common import *
from hl import SymN, F as HF, tt as Htt
import sympy as sp
T = load_table()
PARTS = {n: list(partitions(n)) for n in range(0, 8)}
s_, t_ = ST, TT
def ex(p): return p.as_expr()
def at_t(poly, v): return poly.compose(t, R(v))

# ---------- Thm 2.2 (DS) ----------
main = Tally('Thm 2.2 DS'); neg = Tally('neg DS')
for n in range(1, 8):
    for lam in PARTS[n]:
        for mu in PARTS[n]:
            cc = c(T, lam, mu)
            if not dom(mu, lam):
                main.add(cc == 0, f'support {lam}->{mu}'); continue
            if mu == lam:
                main.add(cc == s**nfun(lam), f'diag {lam}'); continue
            ok = cc != 0 and at_t(cc, 0) != 0 or cc != 0
            ok = ok and cc.compose(s, R(1)) == 0 and val_s(cc) == nfun(mu)
            c0 = at_t(cc, 0)
            ok = ok and c0 != 0 and val_s(c0) == nfun(mu) and min((i for (i, j), _ in c0.terms())) == nfun(mu) \
                 and sum(co for (i, j), co in c0.terms() if i == nfun(mu)) == 1
            main.add(ok, f'{lam}->{mu}')
            neg.add(val_s(cc) == nfun(lam), f'neg {lam}->{mu}')   # perturbed: valuation n(lam) instead of n(mu)
report(main, neg)

# ---------- Remark 2.3 (s->0 edge: d = e->HL transition) ----------
main = Tally('Rem 2.3 d=M(e,P)'); neg = Tally('neg RemH')
main2 = Tally('Rem 2.3 d(1)=#0-1 matrices');
for n in range(1, 8):
    S = SymN(n)
    for lam in PARTS[n]:
        MeP = S.expand_in_P(S.e_to_m({lam: HF(1)}))
        for mu in PARTS[n]:
            d = sum((co * t**j for (i, j), co in c(T, lam, mu).terms() if i == nfun(mu)), R(0))
            mval = MeP.get(conj(mu), HF(0)).as_expr()
            rhs = sp.cancel(t_**(nfun(conj(mu)) - nfun(conj(lam))) * mval.subs(t_, 1 / t_))
            main.add(sp.expand(ex(d) - rhs) == 0, f'{lam},{mu}: {ex(d)} vs {rhs}')
            negr = sp.cancel(mval.subs(t_, 1 / t_))  # perturbed: drop the t-power
            neg.add(sp.expand(ex(d) - negr) == 0)
            main2.add(ex(d).subs(t_, 1) == S.Em[S.idx[lam], S.idx[conj(mu)]], f'd(1) {lam},{mu}')
report(main, neg)
report(main2, neg)

# ---------- Remark 3.6 / proof of Thm 3.4: omega e*_lam|_{t=0} = s^{n(lam)} Q'_lam(x;1/s) ----------
main = Tally('Rem 3.6 t=0 edge (omega e*|t=0 = s^n Q\'(x;1/s))'); neg = Tally('neg t0')
for n in range(1, 8):
    S = SymN(n); K = S.kostka_foulkes()
    for lam in PARTS[n]:
        lhs = {}; rhs = {}; rneg = {}
        for mu, cc in T[lam].items():
            c0 = ex(at_t(cc, 0))
            for j, q in enumerate(S.B):
                if S.Em[S.idx[mu], j]: lhs[q] = lhs.get(q, 0) + c0 * S.Em[S.idx[mu], j]
        for nu in PARTS[n]:
            kf = K[nu].get(lam, HF(0)).as_expr()
            coef = s_**nfun(lam) * kf.subs(t_, 1 / s_)
            coefn = s_**nfun(lam) * kf.subs(t_, s_)   # perturbed: u=s instead of 1/s
            for j, q in enumerate(S.B):
                v = S.Sm[S.idx[conj(nu)], j]
                if v: rhs[q] = rhs.get(q, 0) + coef * v; rneg[q] = rneg.get(q, 0) + coefn * v
        main.add(all(sp.simplify(lhs.get(q, 0) - rhs.get(q, 0)) == 0 for q in S.B), f'{lam}')
        neg.add(all(sp.simplify(lhs.get(q, 0) - rneg.get(q, 0)) == 0 for q in S.B))
report(main, neg)

# ---------- Thm 3.4 block valuation law ----------
main = Tally('Thm 3.4 block law'); neg = Tally('neg blocklaw')
tvals = [QQ(-2), QQ(-1), QQ(2), QQ(1, 2), QQ(1), QQ(3), QQ(-1, 3)]
for n in range(1, 8):
    for lam in PARTS[n]:
        for mu in PARTS[n]:
            cc = c(T, lam, mu); k = kappa(lam, mu); L = len(lam)
            if k == NEG:
                main.add(cc == 0, f'{lam}->{mu} kappa=-inf but c!=0'); continue
            # lower bound for all t (generic + specialisations)
            ok = val_sm1(cc) >= L - k
            for tv in tvals:
                ct = at_t(cc, tv)
                ok = ok and (ct == 0 or val_sm1(ct) >= L - k)
            if mu != lam:
                ok = ok and val_sm1(cc) == L - k
                c0 = at_t(cc, 0)
                ok = ok and c0 != 0 and val_sm1(c0) == L - k
                ld = coeff_sm1(c0, L - k)
                ld = ld.as_expr()
                ok = ok and ld.is_Integer and (-1)**(L - k) * ld >= 1
                neg.add(val_sm1(cc) == L - len(mu))   # perturbed: naive l(lam)-l(mu)
            main.add(ok, f'{lam}->{mu} k={k} v={val_sm1(cc)}')
report(main, neg)

# ---------- Thm 3.7 coarsenings: tight histories x W_k(J) (W from its printed DEFINITION via engine) ----------
E7 = Engine(7, seed=77)
Wcache = {}
def W_def(k, J):
    key = (k, tuple(J))
    if key in Wcache: return Wcache[key]
    p = len(J)
    def Ep(F, deg):  # E_k^{(p)} F for s-free F
        G = E7.E(k, F, deg)
        return {mu: coeff_sm1(v, p) for mu, v in G.items() if coeff_sm1(v, p) != 0}
    def Cn(i, F, deg):  # nested commutator [..[E^{(p)}, e_j1],..,e_ji] applied to F
        if i == 0: return Ep(F, deg)
        j = J[i - 1]
        a = Cn(i - 1, mul_e({(j,): R(1)}, F), deg + j)
        b = mul_e({(j,): R(1)}, Cn(i - 1, F, deg))
        out = dict(a)
        for mu, v in b.items(): out[mu] = out.get(mu, R(0)) - v
        return {m: v for m, v in out.items() if v != 0}
    res = Cn(p, {(): R(1)}, 0)
    n = k + sum(J)
    Wcache[key] = res.get((n,), R(0))
    return Wcache[key]

def W_printed(k, J, p_is='count'):
    n = k + sum(J); p = len(J) if p_is == 'count' else sum(J)
    r = qint(n) * R(1)
    pr = R(1)
    for j in J: pr *= qint(k, t**j)
    num = (-1)**p * qint(n) * pr
    return ex(num) / ex(qint(k))

def histories(lam):
    """tight histories processing lam_ell,...,lam_1; blocks labelled; yields (final sizes, list of (k,J))"""
    order = list(reversed(lam))
    def rec(i, blocks, merges):
        if i == len(order):
            yield tuple(sorted(blocks, reverse=True)), merges; return
        k = order[i]
        yield from rec(i + 1, blocks + [k], merges)
        for p in range(1, len(blocks) + 1):
            for S in itertools.combinations(range(len(blocks)), p):
                J = [blocks[j] for j in S]
                nb = [blocks[j] for j in range(len(blocks)) if j not in S] + [k + sum(J)]
                yield from rec(i + 1, nb, merges + [(k, tuple(sorted(J, reverse=True)))])
    yield from rec(0, [], [])

main = Tally('Thm 3.7 merge-history expansion'); neg = Tally('neg coarse')
mainW = Tally('Thm 3.8 W_k(J) closed form (p=#J) vs definition'); negW = Tally('neg W (p=sum J literal-alt)')
for n in range(2, 8):
    for lam in PARTS[n]:
        for mu in PARTS[n]:
            if mu == lam or not is_coarsening(lam, mu): continue
            m = len(lam) - len(mu); cc = c(T, lam, mu)
            tot = 0; totneg = 0; cnt = 0
            for fin, merges in histories(lam):
                if fin != mu: continue
                cnt += 1
                w = 1; wn = 1
                for (k, J) in merges:
                    wd = ex(W_def(k, J)); w *= wd
                    wn *= wd.subs(t_, 0) if False else wd
                tot += w
            ok = val_sm1(cc) == m and sp.simplify(ex(coeff_sm1(cc, m)) - tot) == 0
            ok = ok and ex(coeff_sm1(at_t(cc, 0), m)) == (-1)**m * cnt and cnt > 0
            main.add(ok, f'{lam}->{mu}')
            neg.add(ex(coeff_sm1(at_t(cc, 0), m)) == cnt)  # perturbed: drop (-1)^m
for (k, J) in list(Wcache.keys()):
    wd = ex(Wcache[(k, J)])
    mainW.add(sp.simplify(wd - W_printed(k, J, 'count')) == 0, f'W_{k}{J}: {sp.factor(wd)}')
    negW.add(sp.simplify(wd - W_printed(k, J, 'sum')) == 0)
# extra W cases beyond those met in histories
for k in range(1, 5):
    for d in range(1, 8 - k):
        for J in PARTS[d]:
            if (k, J) in Wcache: continue
            wd = ex(W_def(k, J))
            mainW.add(sp.simplify(wd - W_printed(k, J, 'count')) == 0, f'W_{k}{J}')
            negW.add(sp.simplify(wd - W_printed(k, J, 'sum')) == 0)
report(main, neg); report(mainW, negW)

# sign consequence of Thm 3.8: (-1)^m Lead > 0 for every real t>0 (coarsenings)
main = Tally('Thm 3.8 sign (-1)^m Lead>0 on t>0'); neg = Tally('neg sign')
for n in range(2, 8):
    for lam in PARTS[n]:
        for mu in PARTS[n]:
            if mu == lam or not is_coarsening(lam, mu): continue
            m = len(lam) - len(mu); Ld = sp.Poly(ex(coeff_sm1(c(T, lam, mu), m)), t_)
            pos_roots = [r for r in sp.real_roots(Ld) if r > 0]
            main.add(not pos_roots and (-1)**m * Ld.eval(1) > 0, f'{lam}->{mu}')
            neg.add(not pos_roots and (-1)**(m + 1) * Ld.eval(1) > 0)
report(main, neg)

# ---------- Thm 4.1 full-merge lead + Cor 4.2 ----------
def conn_sum_brute(ws, L):
    edges = [(i, j) for i in range(L) for j in range(i + 1, L)]
    tot = 0
    for r in range(len(edges) + 1):
        for H in itertools.combinations(edges, r):
            # connected spanning?
            par = list(range(L))
            def f(a):
                while par[a] != a: a = par[a]
                return a
            for (i, j) in H: par[f(i)] = f(j)
            if len({f(i) for i in range(L)}) == 1:
                tot += sp.Mul(*[ws[e] for e in H])
    return sp.expand(tot)
def conn_sum_rec(ws, V):
    # connected part via fixing the min vertex: f(V) = sum_{S ∋ v0} conn(S) f(V\S), f(S)=prod_{e in S}(1+w_e)
    V = tuple(V)
    @lru_cache(None)
    def fS(S): return sp.Mul(*[1 + ws[(i, j)] for i in S for j in S if i < j])
    @lru_cache(None)
    def conn(S):
        if len(S) == 1: return sp.Integer(1)
        v0, rest = S[0], S[1:]; tot = fS(S)
        for r in range(len(rest)):
            for U in itertools.combinations(rest, r):
                A = (v0,) + U; B = tuple(x for x in rest if x not in U)
                tot -= conn(A) * fS(B)
        return sp.expand(tot)
    return conn(V)

main = Tally('Thm 4.1 full-merge lead'); neg = Tally('neg G')
mainK = Tally('Thm 4.1 K = connected part of t^{e2(lam)} (brute graph sum vs cumulant)')
for n in range(2, 8):
    for lam in PARTS[n]:
        L = len(lam)
        if L < 2: continue
        ws = {(i, j): t_**(lam[i] * lam[j]) - 1 for i in range(L) for j in range(i + 1, L)}
        Krec = conn_sum_rec(ws, range(L))
        if L <= 5:
            Kb = conn_sum_brute(ws, L); mainK.add(sp.expand(Kb - Krec) == 0, f'{lam}')
        pref = (1 - t_**n) / sp.Mul(*[1 - t_**l for l in lam])
        lead = ex(coeff_sm1(c(T, lam, (n,)), L - 1))
        main.add(val_sm1(c(T, lam, (n,))) == L - 1 and sp.simplify(lead - pref * Krec) == 0, f'{lam}')
        wsn = {(i, j): t_**(lam[i] + lam[j]) - 1 for i in range(L) for j in range(i + 1, L)}
        neg.add(sp.simplify(lead - pref * conn_sum_rec(wsn, range(L))) == 0)
report(main, neg); report(mainK, neg)

# Cor 4.2
def tree_inversions(n):
    # labelled trees on {1..n} rooted at 1 via Pruefer codes; inv = #{(i,j): i<j, j ancestor of i}
    tot = sp.Integer(0)
    if n == 1: return sp.Integer(1)
    if n == 2: return sp.Integer(1)
    for code in itertools.product(range(1, n + 1), repeat=n - 2):
        deg = [1] * (n + 1)
        for x in code: deg[x] += 1
        edges = []; code = list(code)
        for x in code:
            leaf = min(i for i in range(1, n + 1) if deg[i] == 1)
            edges.append((leaf, x)); deg[leaf] -= 1; deg[x] -= 1
        u, v = [i for i in range(1, n + 1) if deg[i] == 1]; edges.append((u, v))
        adj = {i: [] for i in range(1, n + 1)}
        for a, b in edges: adj[a].append(b); adj[b].append(a)
        par = {1: None}; st = [1]
        while st:
            a = st.pop()
            for b in adj[a]:
                if b not in par: par[b] = a; st.append(b)
        inv = 0
        for i in range(2, n + 1):
            j = par[i]
            while j is not None:
                if j > i: inv += 1
                j = par[j]
        tot += t_**inv
    return sp.expand(tot)
main = Tally('Cor 4.2 (a)-(d)'); neg = Tally('neg CorG')
for n in range(2, 8):
    lead = ex(coeff_sm1(c(T, (1,) * n, (n,)), n - 1))
    In = tree_inversions(n)
    main.add(sp.expand(lead - (-1)**(n - 1) * ex(qint(n)) * In) == 0, f'(a) n={n}')
    neg.add(sp.expand(lead - (-1)**(n - 1) * In) == 0)
    for lam in PARTS[n]:
        L = len(lam)
        if L < 2: continue
        lead = ex(coeff_sm1(c(T, lam, (n,)), L - 1))
        main.add(lead.subs(t_, 0) == (-1)**(L - 1) * sp.factorial(L - 1), f'(b) {lam}')
        neg.add(lead.subs(t_, 0) == (-1)**(L) * sp.factorial(L - 1))
        main.add(sp.limit(lead, t_, 1) == (-1)**(L - 1) * n**(L - 1), f'(c) {lam}')
        neg.add(sp.limit(lead, t_, 1) == (-1)**(L - 1) * n**(L - 2))
        P = sp.Poly(sp.expand(lead), t_)
        main.add(not [r for r in sp.real_roots(P) if r > 0] and sp.sign(P.eval(1)) == (-1)**(L - 1), f'(d) {lam}')
report(main, neg)

# Remark after Cor 4.2: separator J (only the star side is checkable here)
Lf = lambda lam: ex(coeff_sm1(c(T, lam, (sum(lam),)), len(lam) - 1))
J = sp.simplify(Lf((1, 1, 1, 1)) * Lf((2, 2)) / Lf((2, 1, 1))**2)
Jp = (t_**2 + 1) * (t_**3 + 3 * t_**2 + 6 * t_ + 6) / ((t_ + 1) * (t_**2 + t_ + 2)**2)
a_, g4, f1, f2 = sp.symbols('a g4 f1 f2')
gauge = lambda lam: a_**len(lam) * g4 * sp.Mul(*[{1: f1, 2: f2}[x] for x in lam])
Jg = sp.simplify(gauge((1, 1, 1, 1)) * gauge((2, 2)) / gauge((2, 1, 1))**2)
main = Tally('Rem 4.3 separator J_star'); neg = Tally('neg J')
main.add(sp.simplify(J - Jp) == 0, f'J={sp.factor(J)}'); main.add(Jg == 1, 'gauge'); main.add(Jp.subs(t_, 0) == sp.Rational(3, 2)); main.add(Jp.subs(t_, 1) == 1)
neg.add(sp.simplify(J - Jp.subs(t_, -t_)) == 0)
report(main, neg)
# "first gauge invariant appears at n=4": exponent lattice rank test for n<=3 (l>=2 leads)
def expvec(lam, nmax):
    v = [len(lam)] + [1 if sum(lam) == k else 0 for k in range(2, nmax + 1)] + [lam.count(k) for k in range(1, nmax + 1)]
    return v
rows = [expvec(l, 3) for n in (2, 3) for l in PARTS[n] if len(l) >= 2]
print('   gauge: n<=3 leads exponent matrix rank', sp.Matrix(rows).rank(), 'of', len(rows), '(full rank => no invariant at n<=3)')

# ---------- Thm 4.4 (F) and Thm 4.5 (block multiplicativity) ----------
def lead_full(sub):
    sub = srt(sub)
    if len(sub) == 1: return sp.Integer(1)
    return ex(coeff_sm1(c(T, sub, (sum(sub),)), len(sub) - 1))
main = Tally('Thm 4.4 coarsening leads factor'); neg = Tally('neg F')
for n in range(2, 8):
    for lam in PARTS[n]:
        for mu in PARTS[n]:
            if mu == lam or not is_coarsening(lam, mu): continue
            m = len(lam) - len(mu)
            tot = 0; totn = 0
            for pi in set_partitions(list(range(len(lam)))):
                if srt(tuple(sum(lam[i] for i in B) for B in pi)) != mu: continue
                pr = sp.Mul(*[lead_full(tuple(lam[i] for i in B)) for B in pi]); tot += pr
                totn += pr * (2 if len(pi) > 1 and len(pi) < len(lam) else 1)  # perturbed: interleaving factor 2
            ld = ex(coeff_sm1(c(T, lam, mu), m))
            main.add(sp.simplify(ld - tot) == 0, f'{lam}->{mu}')
            neg.add(sp.simplify(ld - totn) == 0)
report(main, neg)

main = Tally('Thm 4.5 block multiplicativity'); neg = Tally('neg blockmult')
mainK1 = Tally('Thm 4.5 side claim: every factor has kappa(lam_C,nu^C)=1')
for n in range(1, 8):
    for lam in PARTS[n]:
        L = len(lam)
        for mu in PARTS[n]:
            k = kappa(lam, mu)
            if k == NEG or k < 1: continue
            m = L - k; tot = 0; totn = 0
            for pi in set_partitions(list(range(L))):
                if len(pi) != k: continue
                blocks = [srt(tuple(lam[i] for i in B)) for B in pi]
                choices = [[nu for nu in PARTS[sum(b)] if dom(nu, b)] for b in blocks]
                for tup in itertools.product(*choices):
                    if srt(sum(tup, ())) != mu: continue
                    pr = 1
                    for b, nu in zip(blocks, tup):
                        pr *= ex(coeff_sm1(c(T, b, nu), len(b) - 1))
                        mainK1.add(kappa(b, nu) == 1, f'{lam}->{mu}: block {b}->{nu} kappa={kappa(b, nu)}')
                    tot += pr
            totn = tot * sp.factorial(k)   # perturbed: ordered blocks (compositions) instead of set partitions
            ld = ex(coeff_sm1(c(T, lam, mu), m))
            main.add(sp.simplify(ld - tot) == 0, f'{lam}->{mu}')
            neg.add(sp.simplify(ld - totn) == 0)
report(main, neg); report(mainK1, Tally('n/a'))

# ---------- Thm 5.1 Box Complement, Cor 5.2 Column Lemma ----------
main = Tally('Thm 5.1 Box Complement'); neg = Tally('neg box')
for L in range(1, 5):
    for N in range(1, 8):
        for d in range(0, 8):
            if N * L - d < 0 or N * L - d > 7: continue
            cand = [p for p in PARTS[d] if len(p) <= L and (not p or p[0] <= N)]
            for lam in cand:
                for mu in cand:
                    lc = srt(tuple(N - x for x in (list(lam) + [0] * (L - len(lam)))[::-1]))
                    mc = srt(tuple(N - x for x in (list(mu) + [0] * (L - len(mu)))[::-1]))
                    lhs = c(T, lc, mc) if lc else (R(1) if not mc else R(0))
                    base = c(T, lam, mu) if lam else (R(1) if not mu else R(0))
                    e = N * L * (L - 1) // 2 - (L - 1) * d
                    # exponent may be negative (Laurent in s): compare lhs*s^{max(0,-e)} with s^{max(0,e)}*base
                    eq = lambda ee: lhs * s**max(0, -ee) == s**max(0, ee) * base
                    main.add(eq(e), f'L={L} N={N} {lam}->{mu} e={e}')
                    if base != 0: neg.add(eq(e + 1))
report(main, neg)
main = Tally('Cor 5.2 Column Lemma (l = max length)'); neg = Tally('neg column')
main2 = Tally('Cor 5.2 Column Lemma with l > max length (zero-padded)')
for n in range(1, 8):
    for lam in PARTS[n]:
        for mu in PARTS[n]:
            L0 = max(len(lam), len(mu))
            for L, tal in ((L0, main), (L0 + 1, main2)):
                if n + L > 7: continue
                l2 = tuple(x + 1 for x in list(lam) + [0] * (L - len(lam)))
                m2 = tuple(x + 1 for x in list(mu) + [0] * (L - len(mu)))
                tal.add(c(T, l2, m2) == s**(L * (L - 1) // 2) * c(T, lam, mu), f'{lam}->{mu} L={L}')
                if tal is main and c(T, lam, mu) != 0: neg.add(c(T, l2, m2) == s**(L * (L + 1) // 2) * c(T, lam, mu))
report(main, neg); report(main2, neg)
ex_lead = ex(qint(3)) * (t_ + 2)
main = Tally('Cor 5.2 example leads [3](t+2)');
for lam, mu in [((2, 2, 2), (3, 3)), ((2, 2, 2), (4, 1, 1)), ((1, 1, 1), (3,))]:
    cc = c(T, lam, mu); main.add(val_sm1(cc) == 2 and sp.expand(ex(coeff_sm1(cc, 2)) - ex_lead) == 0, f'{lam}->{mu}: {sp.factor(ex(coeff_sm1(cc, 2)))}')
report(main, Tally('none'))

# ---------- Cor 5.4 Pieri (all s) + Hikita k=1 rule quoted in Sec. 2 ----------
u_ = sp.Symbol('u')
def pi_a(a, j):
    f = sp.Mul(*[(1 + s_ * u_ * t_**m) / (1 + u_ * t_**m) for m in range(a)])
    return sp.series(f, u_, 0, j + 1).removeO().coeff(u_, j) if j > 0 else sp.Integer(1)
def Cab(a, b):
    if a == 0 or b == 0: return sp.Integer(1)
    return (-1)**b * ex(qint(a + b)) / ex(qint(a)) * pi_a(a, b)
main = Tally('Cor 5.4 Pieri e_k*e_r'); neg = Tally('neg pieri')
for k in range(1, 7):
    for r in range(1, 8 - k):
        lam = srt((k, r)); row = T[lam]
        pred = {}; predn = {}
        for y in range(0, min(k, r) + 1):
            key = srt((k + r - y, y)); v = s_**y * Cab(k - y, r - y)
            pred[key] = pred.get(key, 0) + v; predn[key] = predn.get(key, 0) + Cab(k - y, r - y)
        ok = all(sp.simplify(ex(row.get(mu, R(0))) - pred.get(mu, 0)) == 0 for mu in set(row) | set(pred))
        main.add(ok, f'k={k} r={r}')
        neg.add(all(sp.simplify(ex(row.get(mu, R(0))) - predn.get(mu, 0)) == 0 for mu in set(row) | set(predn)))
report(main, neg)
main = Tally('Sec 2 Hikita Pieri e_1*e_r=(1-s)[r+1]e_{r+1}+s e_1e_r'); neg = Tally('neg hik')
for r in range(1, 7):
    row = T[srt((1, r))]
    main.add(sp.expand(ex(row.get((r + 1,), R(0))) - (1 - s_) * ex(qint(r + 1))) == 0 and row.get(srt((r, 1)), R(0)) == s and len(row) == 2, f'r={r}')
    neg.add(sp.expand(ex(row.get((r + 1,), R(0))) - (1 - s_) * ex(qint(r))) == 0)
report(main, neg)

# ---------- Open problems: computed side-claims ----------
cnt = 0
for lam in PARTS[6]:
    for mu in PARTS[6]:
        if mu == lam or not dom(mu, lam): continue
        k = kappa(lam, mu); ct = at_t(c(T, lam, mu), QQ(-2))
        if ct == 0 or val_sm1(ct) > len(lam) - k: cnt += 1
print(f"RESULT [Open problem 2: #pairs n=6 with val>l-kappa at t=-2] {'PASS' if cnt == 7 else 'FAIL'}: found {cnt} (printed: seven)")
ok = True; vals = []
for n in range(3, 8):
    for lam in PARTS[n]:
        if len(lam) != 3 or kappa(lam, (n - 1, 1)) != 1: continue
        v = ex(coeff_sm1(c(T, lam, (n - 1, 1)), 2)).subs(t_, 1); vals.append((lam, v)); ok = ok and v == 2 * n * n - 6 * n + 3
print(f"RESULT [Open problem 3: Lead(1)=2n^2-6n+3, l=3, mu=(n-1,1), kappa=1] {'PASS' if ok and vals else 'FAIL'}: cases={len(vals)} (n<=7)")
