"""K0 kill test (Wake 249): Psi_{kappa,0}(M) vs Ikeda-Naruse GQ_lambda (1112.5223 Def 2.1 / (2.2),
= Motegi 2610.15344 Def 2.3 (2.17)), beta = -1.
GQ side implemented independently from the PRINTED formula:
  x(+)y = x+y-xy,  x(-)y = (x-y)/(1-y),  [[x|b]]^k = (x(+)x) prod_{j=1}^{k-1} (x(+)b_j)
  GQ_lam(x_1..x_n|b) = 1/(n-r)! sum_{w in S_n} w[ prod_{i<=r} [[x_i|b]]^{lam_i} prod_{i<=r} prod_{j=i+1}^n (x_i(+)x_j)/(x_i(-)x_j) ]
lam strictly DECREASING.  Engine side: prove248/psi.py (T0 definition).
"""
import sys, itertools, random, math
import sympy as sp
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/prove248')
from psi import Psi, W   # engine side only

def op(x, y): return x + y - x*y
def om(x, y): return (x - y)/(1 - y)

def GQ(lam, xs, b):
    n, r = len(xs), len(lam)
    tot = 0
    for p in itertools.permutations(range(n)):
        y = [xs[i] for i in p]
        t = sp.Integer(1)
        for i in range(r):
            t *= op(y[i], y[i])
            for k in range(1, lam[i]):           # k = 1 .. lam_i - 1
                t *= op(y[i], b[k-1])
            for j in range(i+1, n):
                t *= op(y[i], y[j]) / om(y[i], y[j])
        tot += t
    return tot / math.factorial(n - r)

def wv(k): return W[k-1]

# dictionaries: variables x^IN_m for m in M, parameters b_k
DICTS = {
  'X=1/w,B=w (x=1-1/w_m, b=1-w_k)': (lambda m: 1 - 1/wv(m), lambda k: 1 - wv(k)),
  'CONTROL X=w (x=1-w_m, b=1-w_k)': (lambda m: 1 - wv(m),   lambda k: 1 - wv(k)),
  'CONTROL B=1/w (x=1-1/w_m, b=1-1/w_k)': (lambda m: 1 - 1/wv(m), lambda k: 1 - 1/wv(k)),
}
SHIFTS = {'lam=rev(kappa)': 0, 'lam=rev(kappa)-1 (connection file)': 1}

def rand_point(seed):
    rnd = random.Random(seed)
    return {W[i]: sp.Rational(rnd.randint(2, 97), rnd.randint(2, 97)) for i in range(len(W))}

def compare(kap, M, th=None, symbolic=False):
    r = len(kap); th = th or (0,)*r
    eng = Psi(kap, th, M)
    out = {}
    for dn, (xf, bf) in DICTS.items():
        for sn, sh in SHIFTS.items():
            lam = tuple(sorted((k - sh for k in kap), reverse=True))
            if lam[-1] <= 0: out[(dn, sn)] = 'n/a (zero part)'; continue
            xs = [xf(m) for m in M]; b = [bf(k) for k in range(1, max(lam) + 1)]
            if symbolic:
                g = GQ(lam, xs, b)
                d = sp.cancel(sp.together(g - eng))
                if d == 0: out[(dn, sn)] = 'EQUAL (symbolic)'
                else:
                    q = sp.factor(sp.cancel(sp.together(g/eng)))
                    out[(dn, sn)] = 'DIFF; ratio=' + (str(q) if len(str(q)) < 80 else '<non-constant>')
            else:
                oks = []; ratios = []
                for s in range(3):
                    P = rand_point(1000*s + sum(kap) + 7*len(M))
                    g = GQ(lam, [sp.nsimplify(x.subs(P)) for x in xs], [bb.subs(P) for bb in b])
                    e = eng.subs(P)
                    oks.append(sp.simplify(g - e) == 0); ratios.append(sp.nsimplify(g/e) if e != 0 else None)
                out[(dn, sn)] = 'EQUAL (3 rational pts)' if all(oks) else f'DIFF ratios={ratios}'
    return out

if __name__ == '__main__':
    cases = [  # (kappa increasing, M, theta, symbolic?)
        ((1, 2), (1, 2, 3, 4), None, True),
        ((2, 3), (2, 3, 4, 5), None, True),
        ((1, 3), (1, 2, 3, 4), None, True),
        ((2, 4), (2, 3, 4, 5), None, True),     # generic (non-consecutive, no 1)
        ((2, 5), (1, 3, 5, 6), None, False),    # generic, M not an interval
        ((3, 4), (2, 4, 5, 7), None, False),
        ((2,), (1, 2, 4, 5), None, False),      # r=1 sanity
        ((1, 3, 4), (1, 2, 3, 4, 5), None, False),  # r=3
        # twist control: theta != 0 compared against untwisted GQ
        ((2, 4), (2, 3, 4, 5), (1, 0), True),
        ((2, 4), (2, 3, 4, 5), (1, 1), True),
    ]
    for kap, M, th, sym in cases:
        print(f'=== kappa={kap} M={M} theta={th or (0,)*len(kap)} ===', flush=True)
        for k, v in compare(kap, M, th, sym).items():
            print(f'  {k[0]:42s} | {k[1]:36s} | {v}', flush=True)
