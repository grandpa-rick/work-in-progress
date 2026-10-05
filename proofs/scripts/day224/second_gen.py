"""[(s-1)^2] e*_lambda for general l, numeric t, e-basis dicts. E2_a via polarization (order 2): E2_a(prod g) =
sum_i (prod_{j!=i}) E2_a(g_i) + sum_{i<j} (prod others) Gam_a(g_i,g_j). D_a derivation via Thm 7.1."""
from lead2 import *
from functools import lru_cache
import itertools
def emul(*fs):
    out = {(): Fr(1)}
    for f in fs: out = mul(out, f)
    return out
def E2gen(a, F, t, c1, c2):
    out = {}
    for key, v in F.items():
        for i in range(len(key)):
            rest = {key[:i]+key[i+1:]: v}
            out = add(out, mul(rest, c1(a, key[i])))
        for i, j in itertools.combinations(range(len(key)), 2):
            rest = {tuple(k for idx, k in enumerate(key) if idx not in (i, j)): v}
            out = add(out, mul(rest, c2(a, key[i], key[j])))
    return out
def second_general(lam, t):
    lam = list(lam)
    c1 = lru_cache(None)(lambda a, b: toe(E2op(a, b, t), t))
    c2 = lru_cache(None)(lambda a, b, c: toe(Gam(a, b, c, t), t))
    # states: dict order -> e-dict ; process from the right
    cur = {0: {(): Fr(1)}}
    for a in reversed(lam):
        nxt = {}
        for o, F in cur.items():
            # p=0: multiply by e_a
            nxt[o] = add(nxt.get(o, {}), mul(F, {(a,): Fr(1)}))
            if o+1 <= 2: nxt[o+1] = add(nxt.get(o+1, {}), Dder(a, F, t))
            if o+2 <= 2: nxt[o+2] = add(nxt.get(o+2, {}), E2gen(a, F, t, c1, c2))
        cur = nxt
    return cur
