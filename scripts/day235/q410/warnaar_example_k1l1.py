#!/usr/bin/env python3
"""Third route to c_alpha at k=l=1 (M=3), alpha=(1^6): Warnaar §6 printed example
 sum_y t^{3y^2-y}(e'_{n-6y}-e'_{n-6y+2}) = (x1..xn)^{-1} sum_{r,s} t^{s^2} P_{(2^r,1^{2s})}(x;t).
[x1..x6] P_lambda = sum over standard tableaux of prod psi (Macdonald III (5.11'), psi from memory, code
from korff_vs_macdonald.macdonald_psi). Also M=2 check against Warnaar Conj 1.5 form used by Clio Prop 14."""
import sympy as sp
from korff_vs_macdonald import macdonald_psi, t
def std_sum(lam):
    # sum over chains of partitions adding one box at a time, of prod psi_{nu/nu-}
    from functools import lru_cache
    @lru_cache(None)
    def f(nu):
        if sum(nu)==0: return sp.Integer(1)
        tot=0
        for i in range(len(nu)):
            if nu[i]>0 and (i+1==len(nu) or nu[i]>nu[i+1]):
                m=list(nu); m[i]-=1; m=tuple(x for x in m if x>0)
                tot+=f(m+tuple()) * macdonald_psi(tuple(x for x in nu if x>0), m)
        return sp.expand(tot)
    return f(tuple(lam))
n=6; c=0
for s in range(0,4):
    r=(n-2*s)//2
    lam=tuple([2]*r+[1]*(2*s))
    if len(lam)<=n:
        v=std_sum(lam); print(lam, 't^%d'%(s*s), sp.factor(v)); c+=t**(s*s)*v
print("k=l=1, n=6, alpha=(1^6): c =", sp.expand(c), " claim 5-5t^2+t^4:", sp.expand(c-(5-5*t**2+t**4))==0)
