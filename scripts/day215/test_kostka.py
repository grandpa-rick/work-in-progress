"""Compare d_{lam mu}(t) with Kostka-Foulkes product formulas sum_nu A_{nu,lam} B_{nu,mu'} (and conjugate pairings)."""
import sys, itertools
sys.path.insert(0, '/home/agent/projects/scripts/day215')
from analyze import load, conj, parts, dominates, nstat
from kf import KF
def padd(a, b):
    r = [0]*max(len(a), len(b))
    for i, x in enumerate(a): r[i] += x
    for i, x in enumerate(b): r[i] += x
    while len(r) > 1 and r[-1] == 0: r.pop()
    return r
def pmul(a, b):
    r = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): r[i+j] += x*y
    return r
def K(nu, mu, kind):
    if not dominates(nu, mu): return [0]
    p = list(KF(nu, mu)); p += [0]*(nstat(mu)+1-len(p))
    if kind == '1': return [sum(p)]
    if kind == 't': return p
    if kind == 'tilde': return p[::-1]
def norm(p):
    p = list(p)
    while p and p[0] == 0: p.pop(0)
    while len(p) > 1 and p[-1] == 0: p.pop()
    return p
N = int(sys.argv[1]) if __name__ == "__main__" else 0
kinds = ['1', 't', 'tilde']
results = {}
for pair in ['conj', 'same']:
    for ka in kinds:
        for kb in kinds:
            if ka == kb == '1': continue
            exact = shift = tot = 0; first_fail = None
            for n in range(1, N+1):
                D = load(n)
                for lam in parts(n):
                    for mu, d in D[lam].items():
                        d = [int(x) for x in d]; tot += 1
                        s = [0]
                        for nu in parts(n):
                            a = K(conj(nu) if pair == 'conj' else nu, lam, ka); b = K(nu, conj(mu), kb)
                            s = padd(s, pmul(a, b))
                        if s == d: exact += 1
                        elif norm(s) == d: shift += 1
                        elif first_fail is None or (n, len(lam)) < first_fail[0]: first_fail = ((n, len(lam)), lam, mu, d, s)
            print(f'pair={pair} A={ka} (lam side) B={kb} (mu\' side): exact {exact}/{tot}, up-to-shift {shift}; first fail {first_fail[1:] if first_fail else None}')
# lam = 1^n special case: words of content mu' by cocharge
print('lam=1^n, d vs sum_nu f^nu Ktilde_{nu,mu\'}:')
for n in range(1, N+1):
    D = load(n); lam = (1,)*n
    for mu, d in D[lam].items():
        s = [0]
        for nu in parts(n): s = padd(s, pmul(K(nu, lam, '1'), K(nu, conj(mu), 'tilde')))
        print(' ', n, mu, [int(x) for x in d], s, 'MATCH' if s == [int(x) for x in d] else 'diff')
