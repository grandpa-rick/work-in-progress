"""Day 215: checks on d_{lam mu}(t): N[t], d(0)=1, d(1)=#0-1 matrices (rows lam, cols mu'), symmetry, unimodality, degree."""
import pickle, sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0, '/home/agent/projects/scripts/day214')
from ds_pointeval import parts, dominates, nstat
def conj(k): return tuple(sum(1 for p in k if p > i) for i in range(k[0])) if k else ()
def M01(rows, cols):
    if not rows: return int(all(c == 0 for c in cols))
    r, tot = rows[0], 0
    for S in itertools.combinations(range(len(cols)), r):
        if all(cols[j] > 0 for j in S): tot += M01(rows[1:], [c-1 if j in S else c for j, c in enumerate(cols)])
    return tot
def load(n): return pickle.load(open(f'/home/agent/projects/scripts/day215/d_n{n}.pkl', 'rb'))
if __name__ == '__main__':
    N = int(sys.argv[1])
    for n in range(1, N+1):
        D = load(n); npairs = bad_N = bad0 = bad1 = nonsym = nonuni = 0; degfails = []
        for lam in parts(n):
            up = {mu for mu in parts(n) if dominates(mu, lam)}
            assert set(D[lam]) == up, (lam, set(D[lam]) ^ up)
            for mu, p in D[lam].items():
                npairs += 1
                if any(c.denominator != 1 or c < 0 for c in p): bad_N += 1; print('NOT in N[t]', lam, mu, p)
                if p[0] != 1: bad0 += 1; print('d(0)!=1', lam, mu)
                if sum(p) != M01(list(lam), list(conj(mu))): bad1 += 1; print('d(1) != M01', lam, mu, sum(p), M01(list(lam), list(conj(mu))))
                if p != p[::-1]: nonsym += 1
                i = max(range(len(p)), key=lambda j: p[j])
                if not (all(p[j] <= p[j+1] for j in range(i)) and all(p[j] >= p[j+1] for j in range(i, len(p)-1))): nonuni += 1; print('non-unimodal', lam, mu, p)
                deg = len(p)-1; cand = nstat(conj(mu)) - nstat(conj(lam))
                if deg != cand: degfails.append((lam, mu, deg, cand))
        print(f'n={n}: pairs {npairs}, not-in-N[t] {bad_N}, d(0)!=1 {bad0}, d(1)!=M01 {bad1}, non-palindromic {nonsym}, non-unimodal {nonuni}, deg != n(mu\')-n(lam\') : {len(degfails)}', degfails[:5])
