# Independent check of Clio uid438 §3: (2,3) tie via residues. Rick, Wake 252.
# Psi re-implemented from the primitive T^0 definition (d-infinity file §1), one head only, NOT from q532.py.
import sympy as sp, itertools, random
from fractions import Fraction as Fr
z=sp.Symbol('z')
def Psi1(s,t,S,w):
    """Psi_{(s),t}(S) = sum_{m in S, m>=s} (1-x^-2) prod_{k<s}(1-w_k/x) prod_{u in S\\m} (1-1/(w_u x))/(1-w_u/x) * x^t, x=w_m"""
    tot=0
    for m in S:
        if m<s: continue
        x=w[m]; e=(1-x**-2)*x**t
        for k in range(1,s): e*=(1-w[k]/x)
        for u in S:
            if u!=m: e*=(1-1/(w[u]*x))/(1-w[u]/x)
        tot+=e
    return tot
def Phi(g,S,w):
    return sum(g.subs(z,w[m])/sp.Mul(*[w[m]-w[u] for u in S if u!=m]) for m in S)
def P_(S,w): return sp.Mul(*[w[u]*z-1 for u in S])
def D_(S,w): return sp.Mul(*[z-w[u] for u in S])
def eS(S,w): return sp.Mul(*[w[u] for u in S])
def resinf(R): return -sp.residue(R.subs(z,1/z)/z**2, z, 0)
W=sp.symbols('w1:9'); ws={i+1:W[i] for i in range(8)}
def ratpt(n,seed):
    random.seed(seed); return {k:sp.Rational(random.randint(2,97),random.randint(2,97))+k for k in range(1,n+1)}
ok=bad=0
def rep(lab,cond):
    global ok,bad
    if cond: ok+=1
    else: bad+=1
    print(('PASS ' if cond else 'FAIL ')+lab)
# Step 1: one-head collapse Psi_{(s),t}(S) = Phi[z^{t-s} prod_{k<s}(z-w_k) P]/e_S  (symbolic)
for S in [(1,2),(2,3),(1,2,3),(2,3,4),(1,3,4)]:
    for s in (1,2,3):
        for t in (-2,-1,0,1,2):
            lhs=Psi1(s,t,S,ws); rhs=Phi(z**(t-s)*sp.Mul(*[z-ws[k] for k in range(1,s)])*P_(S,ws),S,ws)/eS(S,ws)
            if sp.cancel(lhs-rhs)!=0: rep(f'collapse s={s} t={t} S={S}',False)
rep('step1 collapse symbolic, 5 S x s<=3 x |t|<=2 (75 cases)',True if bad==0 else False)
# Step 2: w1 cancels: tie <=> Phi[P/z] = e - 1/e  (symbolic, algebraic equivalence)
for S in [(1,2),(2,3),(1,2,3),(2,3,4),(1,2,3,4)]:
    e=eS(S,ws)
    tie=1-Psi1(2,0,S,ws)-(e**-2+ws[1]*Psi1(1,-1,S,ws))
    red=(e-1/e-Phi(P_(S,ws)/z,S,ws))/e   # tie*? : 1-Psi20-e^-2-w1Psi1m1 = (e - 1/e - Phi[P/z])/e
    rep(f'step2 tie-expr == (e-1/e-Phi[P/z])/e symbolic S={S}', sp.cancel(tie-red)==0)
    rep(f'tie itself symbolic S={S}', sp.cancel(tie)==0)
# Step 3: residues of R=P/(zD): Res_0 = 1/e, Res_inf = -e, residues at w_m sum to Phi[P/z]
for S in [(2,3),(1,2,3),(2,3,4)]:
    R=P_(S,ws)/(z*D_(S,ws)); e=eS(S,ws)
    r0=sp.residue(R,z,0); ri=resinf(sp.cancel(R)); rw=sum(sp.residue(R,z,ws[m]) for m in S)
    rep(f'Res0=1/e S={S}', sp.cancel(r0-1/e)==0)
    rep(f'Resinf=-e S={S}', sp.cancel(ri+e)==0)
    rep(f'sum Res_wm = Phi[P/z] S={S}', sp.cancel(rw-Phi(P_(S,ws)/z,S,ws))==0)
# Exact rational: tie on larger S (incl. Clio's ten S)
Slist=[(2,3,4),(1,2,3,4),(2,3,4,5),(1,2,3,4,5),(3,4,5),(2,4,5),(2,3),(4,),(1,5,7),(2,3,4,5,6),(1,2,3,4,5,6,7),(1,3,6,8)]
for S in Slist:
    for sd in (11,12):
        w=ratpt(9,sd); e=eS(S,w)
        rep(f'tie exact S={S} seed={sd}', 1-Psi1(2,0,S,w)==e**-2+w[1]*Psi1(1,-1,S,w))
# One-head machine e_S Psi_{(s),t} = -Res0 R - Resinf R, R = z^{t-s} prod_{k<s}(z-w_k) P/D
for S in [(2,3,4),(1,2,4),(3,5,6,7)]:
    w=ratpt(8,21)
    for s in (1,2,3,4):
        for t in (-3,-1,0,2,3):
            R=z**(t-s)*sp.Mul(*[z-w[k] for k in range(1,s)])*P_(S,w)/D_(S,w)
            rhs=-sp.residue(R,z,0)-resinf(sp.cancel(R))
            rep(f'machine S={S} s={s} t={t}', sp.nsimplify(0)+sp.cancel(eS(S,w)*Psi1(s,t,S,w)-rhs)==0)
print('--- NEGATIVE CONTROLS (each must FAIL; counted separately) ---')
nc_fire=nc_tot=0
def ctrl(lab,holds):
    global nc_fire,nc_tot
    nc_tot+=1; nc_fire+= (not holds); print(('fires ' if not holds else 'DOES NOT FIRE ')+lab)
for S in [(2,3,4),(1,2,3,4),(1,5,7)]:
    w=ratpt(9,31); e=eS(S,w); ph=Phi(P_(S,w)/z,S,w)
    ctrl(f'drop Res_0 term: Phi[P/z]=e  S={S}', ph==e)
    ctrl(f'Res_0 sign flipped: Phi[P/z]=e+1/e  S={S}', ph==e+1/e)
    ctrl(f'Res_inf off by one power: Phi[P/z]=e^2-1/e  S={S}', ph==e**2-1/e)
    ctrl(f'tie with e^-1 instead of e^-2 S={S}', 1-Psi1(2,0,S,w)==1/e+w[1]*Psi1(1,-1,S,w))
    ctrl(f'tie with w2 instead of w1 S={S}', 1-Psi1(2,0,S,w)==e**-2+w[2]*Psi1(1,-1,S,w))
    ctrl(f'tie with twist -2 on (1) S={S}', 1-Psi1(2,0,S,w)==e**-2+w[1]*Psi1(1,-2,S,w))
    # machine with a perturbed residue: drop Res_inf
    R=z**(0-2)*(z-w[1])*P_(S,w)/D_(S,w)
    ctrl(f'machine dropping Res_inf S={S}', e*Psi1(2,0,S,w)==-sp.residue(R,z,0))
print(f'\nSUMMARY: {ok} pass, {bad} fail; controls fired {nc_fire}/{nc_tot}')
