"""Q47: LS12 (arXiv 1203.3237) Cor 3.4 eq (3.6) at lambda=omega_m, type C_m, vs our
A_{kappa rho} = psi((-1)^{|theta|} P_{kappa rho}) (PROVE 259 Thm E).  Brute force, m<=3.
Dictionary used: LS thick flag = G/B^-, X^LS_w <-> our X_{w w0}; z = min rep of coset with z.omega = -v_kappa;
our L = LS L^{omega} (w0=-1 in type C), so our L^{-1} = LS L^{-omega}."""
import itertools, sys
from collections import Counter
from fractions import Fraction as F
m = int(sys.argv[1]) if len(sys.argv) > 1 else 2
CTRL = len(sys.argv) > 2  # negative control: drop the w0 twist (z.omega = +v_kappa)
# signed permutations: w = (perm, signs); w(e_i) = signs[i]*e_{perm[i]}
def act(w, v):
    p, s = w; out = [0]*m
    for i in range(m): out[p[i]] += s[i]*v[i]
    return tuple(out)
W = [(p, s) for p in itertools.permutations(range(m)) for s in itertools.product([1,-1], repeat=m)]
def mul(a, b):  # a*b
    return tuple(tuple(x) for x in zip(*[(a[0][b[0][i]], a[1][b[0][i]]*b[1][i]) for i in range(m)]))
def key(w): return (tuple(w[0]), tuple(w[1]))
W = [key(w) for w in W]
mulk = lambda a, b: key(mul(a, b))
pos = []
for i in range(m):
    for j in range(i+1, m):
        for sg in (1, -1):
            v = [0]*m; v[i] = 1; v[j] = -sg; pos.append(tuple(v))  # e_i - e_j, e_i + e_j
    v = [0]*m; v[i] = 2; pos.append(tuple(v))
def ispos(v):
    for x in v:
        if x: return x > 0
def length(w): return sum(1 for a in pos if not ispos(act(w, a)))
def refl(a):  # s_a as signed perm: v -> v - <v,a^vee> a
    n2 = sum(x*x for x in a)
    imgs = []
    for i in range(m):
        e = [0]*m; e[i] = 1
        c = F(2*sum(e[k]*a[k] for k in range(m)), n2)
        im = [e[k] - c*a[k] for k in range(m)]
        j = [k for k in range(m) if im[k] != 0][0]
        imgs.append((j, int(im[j])))
    return key((tuple(x[0] for x in imgs), tuple(x[1] for x in imgs)))
R = [refl(a) for a in pos]
L = {w: length(w) for w in W}
# Bruhat order: transitive closure of w < w t
up = {w: set() for w in W}
for w in W:
    for t in R:
        u = mulk(w, t)
        if L[u] > L[w]: up[w].add(u)
le = {}
for w in W:
    S = {w}; st = [w]
    while st:
        x = st.pop()
        for y in up[x]:
            if y not in S: S.add(y); st.append(y)
    le[w] = S  # elements >= w
def leq(a, b): return b in le[a]
om = tuple([1]*m)
coset = {w: act(w, om) for w in W}
cos = sorted(set(coset.values()))
minrep = {c: min((w for w in W if coset[w] == c), key=lambda w: L[w]) for c in cos}
def cleq(c, d): return leq(minrep[c], minrep[d])
def ccover(c, d):  # c < d cover in W/W_P
    return cleq(c, d) and c != d and L[minrep[d]] == L[minrep[c]] + 1
def bond(c, d):  # d = s_beta c ; |<c, beta^vee>|
    diff = [d[i]-c[i] for i in range(m)]
    nz = [i for i in range(m) if diff[i]]
    if len(nz) == 1: return 1
    return 2
# LS paths of shape omega: straight, or (s1 > s2) joined by a chain of bond-2 covers (a=1/2)
def half_chain(hi, lo):
    if hi == lo: return True
    return any(ccover(mid, hi) and bond(mid, hi) == 2 and half_chain(mid, lo) for mid in cos if cleq(lo, mid))
paths = [((c,), (F(1),)) for c in cos]
for s1 in cos:
    for s2 in cos:
        if s1 != s2 and cleq(s2, s1) and half_chain(s1, s2):
            paths.append(((s1, s2), (F(1, 2), F(1, 2))))
def endpoint(p): return tuple(sum(a*c[i] for c, a in zip(*p)) for i in range(m))
def down(w, c): return max((v for v in W if coset[v] == c and leq(v, w)), key=lambda v: L[v])
def downp(w, p):
    for c in p[0]: w = down(w, c)
    return w
# our side
def v_of(k): return tuple(-1 if (m - i) in k else 1 for i in range(m))  # position i (0-based) = m+1-p  <=> p = m-i
def D(k): return {(a, a+d) for a, p in enumerate(k) for d in range(p)}
def beta(x):
    a, b = x; d = b - a; v = [0]*m
    if d == 0: v[m-1] = 1
    else: v[m-1-d] += 1; v[m-d] -= 1
    return v
strict = [tuple(sorted(c, reverse=True)) for r in range(m+1) for c in itertools.combinations(range(1, m+1), r)]
def ours(k, r):
    th = sorted(D(k) - D(r))
    if not (D(r) <= D(k)): return None
    out = Counter()
    for g in itertools.product('eo', repeat=len(th)):
        G = dict(zip(th, g))
        if any(x[0] == x[1] and G[x] == 'o' for x in th): continue
        ok = all(not (x != y and x[0] <= y[0] and x[1] <= y[1]) or (G[x] == 'e' and G[y] == 'o') for x in th for y in th)
        if not ok: continue
        e = list(map(F, v_of(k)))
        for x in th:
            if G[x] == 'o':
                b = beta(x); e = [e[i] + b[i] for i in range(m)]
        out[tuple(e)] += 1
    return out, (-1)**len(th)
allok = True; nonP = 0
for k in strict:
    z = minrep[v_of(k)] if CTRL else minrep[tuple(-x for x in v_of(k))]
    lsres = {}
    for w in W:
        if not leq(z, w): continue
        terms = Counter()
        for p in paths:
            if cleq(p[0][0], coset[w]) and downp(w, p) == z:
                terms[tuple(-x for x in endpoint(p))] += 1
        if terms:
            if w != minrep[coset[w]]: nonP += 1
            lsres[w] = (terms, (-1)**(L[w]-L[z]))
    # compare
    for r in strict:
        o = ours(k, r)
        w = minrep[v_of(r)] if CTRL else minrep[tuple(-x for x in v_of(r))]
        ls = lsres.get(w)
        oc = o[0] if (o and o[0]) else None
        if (ls is None) != (oc is None) or (ls and (ls[0] != oc or ls[1] != o[1])):
            allok = False; print('MISMATCH', k, r, ls, o)
    if len(lsres) != sum(1 for r in strict if (ours(k, r) or [0])[0]):
        allok = False; print('extra LS terms', k, len(lsres))
print('m=%d  #W=%d  #LSpaths(omega)=%d  nonminimal-w terms=%d  all entries match: %s' % (m, len(W), len(paths), nonP, allok))
# negative control: claim L^{+omega} (dominant formula) instead, compare count of terms at kappa=(m..1)
