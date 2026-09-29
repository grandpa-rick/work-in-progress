from ek_two_col import *
s, t = s_, t_
br = lambda n: sum(t**i for i in range(n)) if n > 0 else 0
for m in (4, 5, 6):
    r, dt = compute(2, 1, 1, m)
    print(m, f'{dt:.2f}s', {l: sp.factor(c) for l, c in r.items()})
# DS consistency: e2*(e1e1) = (C_2 - (1-s)[2] W_2)/s ; C_2 from DS(2,1,1); W_2 from theorem k=2,r=2
def c(n, j):
    if j < 0 or n - j < 0: return 0
    poch = lambda x, N: sp.Mul(*[1 - x*t**i for i in range(N)])
    al = lambda J: 0 if J < 0 else sp.Mul(*[s - t**i for i in range(1, J+1)]) / sp.Mul(*[1 - t*t**i for i in range(J)])
    return poch(s, n-j)/poch(t, n-j)*(al(j) - s*t**(n-j)*al(j-1))
Fn = lambda n, w: sum(c(n, j)*w**j for j in range(n+1))
def key(*parts): return tuple(sorted([p for p in parts if p > 0], reverse=True))
def add(d, k, v): d[k] = sp.simplify(d.get(k, 0) + v)
W = {}
for b in range(3): add(W, key(b, 4-b), s**b*Fn(2-b, t**(2-b)))
rr = 2
C = {}
add(C, key(rr,1,1), s**3); add(C, key(rr,2), s**2*(1-s)*br(2)); add(C, key(rr+1,1), s*(1-s)*(br(rr+1)+t*br(rr)+s)); add(C, key(rr+2), (1-s)**2*br(rr+2)*(br(rr+1)+s))
pred = {}
for l in set(C) | set(W): add(pred, l, (C.get(l, 0) - (1-s)*br(2)*W.get(l, 0))/s)
r, _ = compute(2, 1, 1, 6)
print('DS(2,1,1) consistency:', all(sp.simplify(r.get(l, 0) - pred.get(l, 0)) == 0 for l in set(r) | set(pred)))
print({l: sp.factor(v) for l, v in pred.items()})
