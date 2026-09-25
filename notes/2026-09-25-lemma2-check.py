# Numerical check of MacBeth UID 292 Lemma 2 in a non-abelian (G,+):
# G = S4 (operation "+" := composition p+q = p∘q, apply q first... any fixed convention),
# D = V4 abelian normal, H = G/D ~ S3.  beta(h1,h2) = -s(h1+h2)+s(h1)+s(h2),
# mu_h(x) = -s(h)+x+s(h); claim: beta'-beta = mu_{h2}(theta(h1)) - theta(h1+h2) + theta(h2).
import itertools, random
G = list(itertools.permutations(range(4)))
def op(p, q): return tuple(p[q[i]] for i in range(4))
def inv(p):
    r = [0]*4
    for i, v in enumerate(p): r[v] = i
    return tuple(r)
e = tuple(range(4))
V4 = [e, (1,0,3,2), (2,3,0,1), (3,2,1,0)]
cos = {}
for g in G:
    key = frozenset(op(g, d) for d in V4); cos.setdefault(key, g)
H = list(cos.keys()); H0 = frozenset(V4)
def hop(a, b): return frozenset(op(x, y) for x in a for y in b)  # coset product
random.seed(1); bad = 0; trials = 0
for _ in range(300):
    s = {h: (e if h == H0 else random.choice(sorted(h))) for h in H}
    sp = {h: (e if h == H0 else random.choice(sorted(h))) for h in H}
    th = {h: op(inv(s[h]), sp[h]) for h in H}
    assert all(t in V4 for t in th.values())
    for h1 in H:
        for h2 in H:
            b = op(op(inv(s[hop(h1,h2)]), s[h1]), s[h2])
            bp = op(op(inv(sp[hop(h1,h2)]), sp[h1]), sp[h2])
            mu = op(op(inv(s[h2]), th[h1]), s[h2])
            rhs = op(op(mu, inv(th[hop(h1,h2)])), th[h2])   # D abelian: order irrelevant
            lhs = op(inv(b), bp)                             # beta' - beta
            trials += 1; bad += (lhs != rhs)
print("Lemma 2 checks:", trials, "failures:", bad)
