"""Independent re-verification of Clio's Theorem B (UID 290) -- Rick, 2026-09-25.
Operators on QH*(Gr_{k,n}) basis = k-subsets S of Z/nZ (bead sites).
Coefficients are exact polynomials in q, stored as {q_exponent: int}.
"""
import itertools, sys
from collections import defaultdict

def padd(a, b, s=1):
    c = dict(a)
    for e, v in b.items():
        c[e] = c.get(e, 0) + s * v
        if c[e] == 0: del c[e]
    return c
def pmul(a, b):
    c = defaultdict(int)
    for e1, v1 in a.items():
        for e2, v2 in b.items():
            c[e1 + e2] += v1 * v2
    return {e: v for e, v in c.items() if v}

class Op:
    """sparse operator: M[S] = {T: poly}; meaning  op(sigma_S) = sum_T poly * sigma_T"""
    def __init__(self, states, M=None):
        self.states = states; self.M = M if M is not None else {S: {} for S in states}
    @staticmethod
    def identity(states, c=None):
        c = c or {0: 1}
        return Op(states, {S: {S: dict(c)} for S in states})
    def __add__(self, o): return self.lin(o, 1)
    def __sub__(self, o): return self.lin(o, -1)
    def lin(self, o, s):
        M = {}
        for S in self.states:
            d = {T: dict(p) for T, p in self.M[S].items()}
            for T, p in o.M[S].items():
                d[T] = padd(d.get(T, {}), p, s)
                if not d[T]: del d[T]
            M[S] = d
        return Op(self.states, M)
    def scale(self, c):  # c: poly
        return Op(self.states, {S: {T: pmul(p, c) for T, p in self.M[S].items() if pmul(p, c)} for S in self.states})
    def after(self, o):
        """self o o  (apply o first, then self)"""
        M = {}
        for S in self.states:
            d = {}
            for U, p in o.M[S].items():
                for T, r in self.M[U].items():
                    d[T] = padd(d.get(T, {}), pmul(p, r))
                    if not d[T]: del d[T]
            M[S] = d
        return Op(self.states, M)
    def __eq__(self, o): return all(self.M[S] == o.M[S] for S in self.states)
    def is_zero(self): return all(not self.M[S] for S in self.states)

def setup(n, k):
    return [frozenset(c) for c in itertools.combinations(range(n), k)]

def gen_a(n, states, i):
    """Postnikov generator a_i: bead i -> i+1 mod n, factor q if i = n-1"""
    M = {}
    for S in states:
        j = (i + 1) % n
        if i in S and j not in S:
            M[S] = {(S - {i}) | {j}: {1 if i == n - 1 else 0: 1}}
        else:
            M[S] = {}
    return Op(states, M)

def word(n, states, A, order):
    """apply generators in the given application order (first element applied first)"""
    op = Op.identity(states)
    for i in order:
        op = A[i].after(op)
    return op

def runs(I, n):
    """maximal cyclic runs of a proper subset I of Z/n, each as list [i, i+1, ...]"""
    I = set(I); out = []
    for i in I:
        if (i - 1) % n not in I:
            r = [i]
            while (r[-1] + 1) % n in I: r.append((r[-1] + 1) % n)
            out.append(r)
    return out

def heap_monomial(n, states, A, run_list, eps_list):
    """a_{I,eps}: for each run with orientation vector eps (eps_s=0: a_s before a_{s+1}),
    build a linear extension; runs commute."""
    order = []
    for r, eps in zip(run_list, eps_list):
        # linear extension of the path poset: generate by topological sort
        m = len(r); before = {x: set() for x in range(m)}
        for s in range(m - 1):
            if eps[s] == 0: before[s + 1].add(s)
            else: before[s].add(s + 1)
        done = []; rem = set(range(m))
        while rem:
            x = min(y for y in rem if before[y] <= set(done))
            done.append(x); rem.remove(x)
        order += [r[x] for x in done]
    return word(n, states, A, order)

def h_op(n, states, A, r):
    """Postnikov h_r = sum over proper I, |I|=r, runs ordered a_s before a_{s+1} (eps=0)"""
    tot = Op(states)
    if r == 0: return Op.identity(states)
    for I in itertools.combinations(range(n), r):
        if len(I) == n: continue
        rl = runs(I, n)
        tot = tot + heap_monomial(n, states, A, rl, [[0] * (len(x) - 1) for x in rl])
    return tot

def e_op(n, states, A, r):
    tot = Op(states)
    if r == 0: return Op.identity(states)
    for I in itertools.combinations(range(n), r):
        if len(I) == n: continue
        rl = runs(I, n)
        tot = tot + heap_monomial(n, states, A, rl, [[1] * (len(x) - 1) for x in rl])
    return tot

def p_heap(n, states, A, e):
    """p_e = R_e(-1) via Clio's orientation sum (5): sum_i sum_eps (-1)^{|eps|} a_{[i,i+e-1],eps}"""
    tot = Op(states)
    for i in range(n):
        r = [(i + s) % n for s in range(e)]
        for eps in itertools.product([0, 1], repeat=e - 1):
            tot = tot + heap_monomial(n, states, A, [r], [list(eps)]).scale({0: (-1) ** sum(eps)})
    return tot

def p_bead(n, states, e):
    """Clio's bead formula (3) at t=-1 (e <= n-1)"""
    M = {}
    for S in states:
        d = {}
        for j in S:
            tgt = (j + e) % n
            if tgt in S: continue
            jumped = sum(1 for s in range(1, e) if (j + s) % n in S)
            qq = 1 if any((j + s) % n == n - 1 for s in range(e)) else 0
            T = (S - {j}) | {tgt}
            d[T] = padd(d.get(T, {}), {qq: (-1) ** jumped})
            if not d[T]: del d[T]
        M[S] = d
    return Op(states, M)

def p_rimhook(n, k, states, e):
    """Rick's formula: alternant a_{beta + e eps_j}, reduced by y^n = (-1)^{k-1} q (any e >= 1)."""
    M = {}
    for S in states:
        beta = sorted(S, reverse=True)
        d = {}
        for idx in range(k):
            g = list(beta); g[idx] += e
            w, g[idx] = divmod(g[idx], n)
            if len(set(g)) < k: continue
            # sign of sorting g into decreasing order
            inv = sum(1 for a in range(k) for b in range(a + 1, k) if g[a] < g[b])
            coeff = (-1) ** inv * ((-1) ** (k - 1)) ** w
            T = frozenset(g)
            d[T] = padd(d.get(T, {}), {w: coeff})
            if not d[T]: del d[T]
        M[S] = d
    return Op(states, M)

def run(n, k, verbose=True):
    states = setup(n, k)
    A = [gen_a(n, states, i) for i in range(n)]
    H = [h_op(n, states, A, r) for r in range(n)]            # h_0..h_{n-1}
    E = [e_op(n, states, A, r) for r in range(n)]            # e_0..e_{n-1}
    P = [None] + [p_heap(n, states, A, e) for e in range(1, n)]
    res = {}
    res['heap==bead'] = all(P[e] == p_bead(n, states, e) for e in range(1, n))
    res['heap==rimhook'] = all(P[e] == p_rimhook(n, k, states, e) for e in range(1, n))
    # Newton in e's with Postnikov's e-operators (independent of h and of R_e):
    PE = [None]
    for m in range(1, n + 1):
        acc = Op(states)
        for i in range(1, m):
            if i < n: acc = acc + E[i].after(PE[m - i]).scale({0: (-1) ** (i - 1)})
        if m < n: acc = acc + E[m].scale({0: (-1) ** (m - 1) * m})
        PE.append(acc)
    res['R_e(-1)==p_e(e-ops)'] = all(P[e] == PE[e] for e in range(1, n))
    res['p_n(e-ops)==(-1)^{k-1}kq'] = PE[n] == Op.identity(states, {1: (-1) ** (k - 1) * k})
    res['rimhook e=n == (-1)^{k-1}kq'] = p_rimhook(n, k, states, n) == Op.identity(states, {1: (-1) ** (k - 1) * k})
    # h_r = 0 as operators for n-k < r <= n-1
    res['h_r=0 for r>n-k'] = all(H[r].is_zero() for r in range(n - k + 1, n))
    # Theorem B
    okB = True
    for r in range(1, n):
        rhs = Op(states)
        for e in range(1, r + 1): rhs = rhs + P[e].after(H[r - e])
        if not (H[r].scale({0: r}) == rhs): okB = False
    res['TheoremB r<=n-1'] = okB
    D = Op(states)
    for e in range(1, n): D = D + P[e].after(H[n - e])
    res['defect==(-1)^{k-1}(n-k)q'] = D == Op.identity(states, {1: (-1) ** (k - 1) * (n - k)})
    # commutation
    res['[p_e,h_j]=0'] = all(P[e].after(H[j]) == H[j].after(P[e]) for e in range(1, n) for j in range(1, n))
    # completed Newton at r=n with h_n := (-1)^{k-1} q, p_n := (-1)^{k-1} k q
    res['completed r=n'] = (D + Op.identity(states, {1: (-1) ** (k - 1) * k})) == Op.identity(states, {1: (-1) ** (k - 1) * n})
    return res

if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    allok = True
    for n in range(3, N + 1):
        for k in range(1, n):
            r = run(n, k)
            bad = [key for key, v in r.items() if not v]
            allok &= not bad
            print(n, k, 'ALL OK' if not bad else 'FAIL: ' + ', '.join(bad), flush=True)
    print('GLOBAL', 'OK' if allok else 'FAILURES')
