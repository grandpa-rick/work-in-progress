"""Is phiinf_I = (vertex monomial) * (phi0_J transported by a signed permutation sigma of W(C_N)) for some J?
sigma: w_i -> w_{pi(i)}^{eps_i}; it acts on fixed points (sign vectors; M = {m: sign_m=+1}) compatibly.
We test BOTH conventions for the point action (M -> sigma M and M -> sigma^{-1} M) and report columns I
for which some (sigma, J) gives a monomial ratio at every M (with matching supports)."""
import sys, pickle, itertools, sympy as sp
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake246')
from common import *
N = int(sys.argv[1]); S = subsets(N); w = sp.symbols('w1:%d' % (N+1)); dim = N*(N+1)//2
A = load(N)
Pinf, _ = pickle.load(open('phis_N%d.pkl' % N, 'rb'))
P0 = {(M, I): A[M, I].subs(q, 0) for M in S for I in S}
def ismono(e):
    e = sp.cancel(e); n, d = sp.fraction(e)
    return e != 0 and len(sp.Poly(n, *w).terms()) == 1 and len(sp.Poly(d, *w).terms()) == 1
def pt(M, pi, eps):   # signs s_i (+1 iff i in M); new sign vector at position pi(i) is eps_i*s_i
    s = {i: (1 if i in M else -1) for i in range(1, N+1)}
    return tuple(sorted(pi[i-1] for i in range(1, N+1) if eps[i-1]*s[i] == 1))
hits = []
for pi in itertools.permutations(range(1, N+1)):
    for eps in itertools.product((1, -1), repeat=N):
        sub = {w[i]: w[pi[i]-1]**eps[i] for i in range(N)}
        for I in S:
            for J in S:
                good = True
                for M in S:
                    a = Pinf[M, I]; b = P0[pt(M, pi, eps), J]
                    b = b.subs(sub, simultaneous=True) if b != 0 else 0
                    if (a == 0) != (b == 0) or (a != 0 and not ismono(a / b)): good = False; break
                if good: hits.append((pi, eps, I, J))
cols = sorted(set(h[2] for h in hits))
print('N=%d: %d (sigma,I,J) matches; columns I matched: %s of %d' % (N, len(hits), cols, len(S)))
for h in hits[:40]: print('  pi=%s eps=%s  phiinf_%s ~ phi0_%s' % h)
