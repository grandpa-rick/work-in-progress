"""Exact rational-function field Q(w_1..w_N) on python-flint fmpq_mpoly, with gcd-reduced fractions.
Drop-in replacement for the FLD of prove247/sy.py (same series/F/coord code)."""
import flint
class Ctx:
    def __init__(s,names):
        s.ctx=flint.fmpq_mpoly_ctx.get(tuple(names),'degrevlex'); s.gens=s.ctx.gens()
        s.one_p=s.ctx.from_dict({(0,)*len(names):1}); s.zero_p=s.one_p-s.one_p
class Fr:
    __slots__=('n','d')
    C=None
    def __init__(s,n,d=None,red=True):
        if d is None: d=Fr.C.one_p
        if red:
            if n==0: d=Fr.C.one_p
            else:
                g=n.gcd(d)
                if g!=1: n=n/g; d=d/g   # exact division
            lc=d.leading_coefficient()
            if lc!=1: n=n/lc; d=d/lc
        s.n=n; s.d=d
    def __eq__(s,o):
        if isinstance(o,Fr): return s.n*o.d==o.n*s.d
        return s.n==o*s.d
    def __hash__(s): return 0
    def __add__(s,o):
        if not isinstance(o,Fr): o=Fr(Fr.C.one_p*o)
        if s.d==o.d: return Fr(s.n+o.n,s.d)
        return Fr(s.n*o.d+o.n*s.d,s.d*o.d)
    __radd__=__add__
    def __neg__(s): return Fr(-s.n,s.d,False)
    def __sub__(s,o): return s+(-o if isinstance(o,Fr) else Fr(Fr.C.one_p*(-o)))
    def __rsub__(s,o): return (-s)+o
    def __mul__(s,o):
        if not isinstance(o,Fr): return Fr(s.n*o,s.d)
        return Fr(s.n*o.n,s.d*o.d)
    __rmul__=__mul__
    def inv(s):
        assert s.n!=0
        return Fr(s.d,s.n)
    def __truediv__(s,o):
        if not isinstance(o,Fr): return Fr(s.n,s.d*o)
        return s*o.inv()
    def __rtruediv__(s,o): return Fr(Fr.C.one_p*o)/s
