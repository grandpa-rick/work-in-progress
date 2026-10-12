# Day 237 (b): test DFK19 Rem 4.8/(4.14) literally: Sigma: X -> X/(t-1), Sigma^{-1}: X -> X(t-1);
# claim (Sigma^{-1} nabla Sigma)|_{t->1/t} P_lam(x;q,t) = t^{-n(lam)} q^{n(lam')} P_lam(x;q,t)  [DFK19 Rem 4.8]
# and the eps-corrected version Sigma_eps: X -> eps X/(t-1)  (eps: p_r -> (-1)^r p_r) -- which is a degree sign,
# and the corrected X -> -X/(t-1) = X/(1-t).
import sys; sys.argv=['x','3']
exec(open('rem63_bridge_check.py').read().split('# ---------------- run')[0])
NMAX=3
for n in range(2, NMAX+1):
    parts=list(partitions(n)); mp=m_in_p(n); S={l:schur(l,n) for l in parts}
    Pq=macdonaldP(n,q,t,mp); Ht={l:gh_Htilde(l,n,S) for l in parts}
    ev={l:t**nn(l)*q**nn(conj(l)) for l in parts}
    def nab(F): a=expand_in(F,Ht,n); return sp.expand(sum(a[l]*ev[l]*Ht[l] for l in parts))
    sig=lambda G: sp.expand(sp.together(G.subs(t,1/t)))
    for name,Sg,Sgi in [('DFK literal X/(t-1)',lambda r:1/(t**r-1),lambda r:(t**r-1)),
                        ('sign-fixed X/(1-t)',lambda r:1/(1-t**r),lambda r:(1-t**r)),
                        ('eps-twist eps X/(t-1)',lambda r:(-1)**r/(t**r-1),lambda r:(-1)**r*(t**r-1))]:
        for l in parts:
            G=sig(Pq[l])                      # (.)|_{t->1/t}: conjugate the operator by t->1/t
            out=sig(pleth(nab(pleth(G,Sg)),Sgi))
            ok=iszero(out - t**(-nn(l))*q**nn(conj(l))*Pq[l], n)
            print(n,name,l,'claim holds' if ok else 'FAILS', flush=True)
