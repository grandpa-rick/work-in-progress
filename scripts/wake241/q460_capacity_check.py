"""Independent check of Clio Q460 (Thm 12) from GGW 2412.18085v1 printed eq:IK-delta1 + eq:delta-sym.
Cut matrix M[P,Q] = prod_{i<=a<j} Delta_{P_i,Q_j}(x_i,x_j) on S=(1^{4a+2}), prefix cols E (even), suffix O.
Rank + projective row count mod p at random x, for several q (generic / special) and weight controls."""
import itertools, random, sys
from sympy import nextprime

def prime_1mod(m, start=2**31):
    p = start
    while True:
        p = nextprime(p)
        if (p - 1) % m == 0: return p
P = prime_1mod(12)   # contains 4th, 3rd, 6th roots of unity
inv = lambda a: pow(a % P, P - 2, P)

def root_of_unity(order):
    for g in range(2, 10**6):
        r = pow(g, (P - 1) // order, P)
        if all(pow(r, order // d, P) != 1 for d in (2, 3) if order % d == 0) and r != 1: return r

# printed eq:IK-delta1 tie-free arc/arc cases (all rows are arcs since every column sum is 1)
def dA(x, y, q): return (y-q*q*x)*(y-q**3*x)*inv((y-x)*(y-q*x))
def dB(x, y, q): return (y-q*q*x)*(q*q*y-x)*inv((y-x)**2)
def dC(x, y, q): return (y-q*q*x)**2*(x-q*y)*inv((y-x)**2*(y-q*x))

def tau(u, v):
    (u1, u2), (v1, v2) = u, v
    assert len({u1, u2, v1, v2}) == 4, "collision -- impossible on S=(1^M)"
    if u2 < v1: return 'A'
    if v2 < u1: return 'AT'
    if u1 < v1 < v2 < u2: return 'B'
    if v1 < u1 < u2 < v2: return 'BT'
    if u1 < v1 < u2 < v2: return 'C'
    if v1 < u1 < v2 < u2: return 'CT'

def Delta(t, x, y, q, mode):
    if mode == 'const':          # control: x-independent letters, distinct per type
        return {'A': 2, 'AT': 3, 'B': 5, 'BT': 5, 'C': 7, 'CT': 11}[t]
    if mode == 'collapseBC':     # control: B,C letters replaced by A (kills trichotomy)
        t = {'B': 'A', 'BT': 'AT', 'C': 'A', 'CT': 'AT'}.get(t, t)
    f = {'A': dA, 'B': dB, 'C': dC}[t.rstrip('T')]
    return f(y, x, q) % P if t.endswith('T') else f(x, y, q) % P   # eq:delta-sym

def ordered_matchings(cols):
    cols = sorted(cols)
    def pm(s):
        if not s: yield []; return
        f = s[0]
        for k in range(1, len(s)):
            for rest in pm(s[1:k] + s[k+1:]): yield [(f, s[k])] + rest
    for m in pm(cols):
        for perm in itertools.permutations(m): yield list(perm)

def rank_mod(rows):
    A = [r[:] for r in rows]; rk = 0; ncol = len(A[0])
    for c in range(ncol):
        piv = next((i for i in range(rk, len(A)) if A[i][c]), None)
        if piv is None: continue
        A[rk], A[piv] = A[piv], A[rk]; iv = inv(A[rk][c])
        A[rk] = [v*iv % P for v in A[rk]]
        for i in range(len(A)):
            if i != rk and A[i][c]:
                f = A[i][c]; A[i] = [(vi - f*vr) % P for vi, vr in zip(A[i], A[rk])]
        rk += 1
    return rk

def proj_count(rows):
    seen = set()
    for r in rows:
        k = next(v for v in r if v); iv = inv(k); seen.add(tuple(v*iv % P for v in r))
    return len(seen)

def run(a, q, mode='IK', split='interleave', seed=1):
    random.seed(seed)
    N, M = 2*a+1, 4*a+2
    if split == 'interleave':
        E = list(range(2, 4*a+1, 2)); O = [c for c in range(1, M+1) if c not in E]
    else:
        E = list(range(1, 2*a+1)); O = list(range(2*a+1, M+1))
    x = [random.randrange(1, P) for _ in range(N)]
    Ps = list(ordered_matchings(E)); Qs = list(ordered_matchings(O))
    rows = []
    for Pm in Ps:
        row = []
        for Qm in Qs:
            v = 1
            for i in range(a):
                for jj in range(a+1):
                    v = v * Delta(tau(Pm[i], Qm[jj]), x[i], x[a+jj], q, mode) % P
            row.append(v)
        rows.append(row)
    return len(Ps), len(Qs), rank_mod(rows), proj_count(rows)

if __name__ == '__main__':
    random.seed(0); qg = random.randrange(2, P)
    qs = {'generic': qg, 'q=-1': P-1, 'q prim 3rd': root_of_unity(3), 'q prim 4th': root_of_unity(4),
          'q prim 6th (q^3=-1)': root_of_unity(6)}
    amax = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    print('p =', P)
    for a in range(1, amax+1):
        for name, q in qs.items():
            print(f'a={a} IK {name:22s} (#P,#Q,rank,proj) =', run(a, q))
        print(f'a={a} control const-letters       ', run(a, qg, mode='const'))
        print(f'a={a} control collapse B,C->A     ', run(a, qg, mode='collapseBC'))
        print(f'a={a} control block split         ', run(a, qg, split='block'))
