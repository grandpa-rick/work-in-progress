from hl import *
def ep(k): return e_p(k)
def Dop(a, b, t):  # D_a(e_b)
    out = {}
    for r in range(1, b+1):
        out = add(out, mul(ep(b-r), T(a, {(r,): Fr(1)}, t)), (-1)**(r-1))
    return out
def c2(r):  # [u^r] e_2(y)
    f = {}
    for i in range(1, r):
        f = add(f, mul({(i,): Fr(1)}, {(r-i,): Fr(1)}), Fr((-1)**r, 2))
    f = add(f, {(r,): Fr(1)}, Fr(-(-1)**r*(r-1), 2))
    return f
def E2op(a, b, t):  # E_a^{(2)}(e_b)
    out = {}
    for r in range(2, b+1):
        out = add(out, mul(ep(b-r), T(a, c2(r), t)))
    return out
def Gam(a, b, c, t):  # bilinear part of E_a^{(2)} on (e_b, e_c)
    out = {}
    for r in range(1, b+1):
        for q in range(1, c+1):
            out = add(out, mul(mul(ep(b-r), ep(c-q)), T(a, {tuple(sorted((r, q), reverse=True)): Fr(1)}, t)), (-1)**(r+q))
    return out
def E2prod(a, b, c, t):
    return add(add(Gam(a, b, c, t), mul(ep(b), E2op(a, c, t))), mul(ep(c), E2op(a, b, t)))
def toe(f, t):
    if not f: return {}
    n = sum(next(iter(f))); return S(n, t).to_e(f)
def L(a, b, t): return (1-t**(a+b))*(t**(a*b)-1)/((1-t**a)*(1-t**b))
def Mthm(k, r, t):  # Thm 7.1 B(e_k,e_r), any order
    if k < r: k, r = r, k
    out = {tuple(sorted((k, r), reverse=True)): Fr(r)} if r > 0 else {}
    for j in range(1, r+1):
        key = tuple(x for x in sorted((k+j, r-j), reverse=True) if x)
        out[key] = out.get(key, 0) + L(k-r+j, j, t)
    return out
