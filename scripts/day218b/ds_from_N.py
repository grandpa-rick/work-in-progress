"""Day 218b: DS (support, lead, val_s c = n(mu), lowest coeff) from Theorem (N).
(N): e*_lam = t^{-n(lam')} N(e_lam), N P_nu = s^{n(nu')} t^{n(nu)} P_nu, P_nu = P_nu(x; q=s, T=1/t).
Checks, symbolic in s,t, all lam |- n <= N:
 (0) e*_lam from (N) == e*_lam from Hikita operator matrices (mats_N5.pkl)  [re-check of N]
 (1) a_{lam nu} := [P_nu] e_lam : support nu <= lam', a_{lam lam'}=1, s-integral (val_s >= 0)
 (2) b_{nu mu} := [e_mu] P_nu : support mu >= nu', b_{nu nu'}=1, s-integral
 (3) val_s c_{lam mu} = n(mu) and [s^{n(mu)}] c = t^{n(mu')-n(lam')} a^{HL}_{lam mu'}(1/t)
 (4) a^{HL}_{lam nu}(tau=0) = K_{nu', lam} > 0 for nu <= lam'  (Schur expansion of e_lam)
 (5) s=1: c_{lam mu}|_{s=1} = delta  ;  P_nu regular at s=1
"""
import sys, pickle, itertools, sympy as sp
sys.path.insert(0,'/home/agent/projects/scripts/day214')
from ek_subset_engine import transpose, parts, dominates
from sympy.utilities.iterables import partitions
s,t,q,T=sp.symbols('s t q T')
N=int(sys.argv[1])
mats=pickle.load(open('/home/agent/projects/scripts/day216/mats_N5.pkl','rb'))
def nn(m): return sum(i*x for i,x in enumerate(m))
def apply(k,vec):
    res={}
    for mu,c in vec.items():
        for nu,d in mats[(k,sum(mu))][mu].items(): res[nu]=res.get(nu,0)+c*d
    return {a:sp.cancel(b) for a,b in res.items() if sp.cancel(b)!=0}
def plist(n): return sorted([tuple(sorted(sum(([k]*v for k,v in p.items()),[]),reverse=True)) for p in partitions(n)])
def m_of(expr,xs):
    P=sp.Poly(sp.expand(expr),*xs); d={}
    for mon,c in P.terms():
        if list(mon)==sorted(mon,reverse=True): d[tuple(a for a in mon if a>0)]=c
    return d
def zee(l):
    from collections import Counter
    r=1
    for k,v in Counter(l).items(): r*=k**v*sp.factorial(v)
    return r
def val(x,var):
    x=sp.cancel(x)
    if x==0: return sp.oo
    nu,de=sp.fraction(x)
    def v(p):
        p=sp.Poly(p,var); return min(m[0] for m in p.monoms())
    return v(nu)-v(de)
def lowcoef(x,var,k):
    return sp.cancel(sp.limit(sp.cancel(x*var**(-k)),var,0)) if False else sp.cancel((sp.cancel(x*var**(-k))).subs(var,0))
bad=0; ok=0
def chk(cond,msg):
    global bad,ok
    if cond: ok+=1
    else: bad+=1; print('FAIL',msg,flush=True)
for n in range(1,N+1):
    xs=sp.symbols(f'x1:{n+1}'); Ps=plist(n)
    pm={l:m_of(sp.prod([sum(x**r for x in xs) for r in l]),xs) for l in Ps}
    Mp=sp.Matrix([[pm[l].get(m,0) for m in Ps] for l in Ps]); Minv=Mp.inv()
    em={l:m_of(sp.prod([sum(sp.prod(c) for c in itertools.combinations(xs,r)) for r in l]),xs) for l in Ps}
    Me=sp.Matrix([[em[l].get(m,0) for m in Ps] for l in Ps]); Meinv=Me.inv()  # m_k = sum Meinv[k,l] e_l
    def ip(f,g):
        fp=[sum(f.get(m,0)*Minv[j,i] for j,m in enumerate(Ps)) for i in range(len(Ps))]
        gp=[sum(g.get(m,0)*Minv[j,i] for j,m in enumerate(Ps)) for i in range(len(Ps))]
        return sum(fp[i]*gp[i]*zee(l)*sp.prod([(1-q**r)/(1-T**r) for r in l]) for i,l in enumerate(Ps))
    P={}
    for lam in Ps:  # lex increasing
        f={lam:sp.Integer(1)}
        for mu in Ps:
            if mu==lam: break
            if dominates(lam,mu):
                c=sp.cancel(ip({lam:1},P[mu])/ip(P[mu],P[mu]))
                for k,v in P[mu].items(): f[k]=sp.cancel(f.get(k,0)-c*v)
        P[lam]={k:sp.factor(v) for k,v in f.items() if v!=0}
    sub={q:s,T:1/t}
    Pm=sp.Matrix([[sp.cancel(sp.sympify(P[nu].get(k,0)).subs(sub)) for k in Ps] for nu in Ps])  # rows P_nu, cols m_k
    B=(Pm*Meinv).applyfunc(sp.cancel)   # B[nu,mu] = [e_mu] P_nu
    A=B.inv().applyfunc(sp.cancel)      # e_lam = sum_nu A[lam,nu] P_nu
    idx={l:i for i,l in enumerate(Ps)}
    # (1),(2)
    for lam in Ps:
        for nu in Ps:
            a=A[idx[lam],idx[nu]]; b=B[idx[nu],idx[lam]]
            if a!=0: chk(dominates(transpose(lam),nu),f'a supp {lam},{nu}'); chk(val(a,s)>=0,f'a s-int {lam},{nu}')
            if b!=0: chk(dominates(lam,transpose(nu)),f'b supp {nu},{lam}'); chk(val(b,s)>=0,f'b s-int {nu},{lam}')
        chk(A[idx[lam],idx[transpose(lam)]]==1 and B[idx[transpose(lam)],idx[lam]]==1,'unitri')
    # (2b) regular at s=1 & P_nu(1,tau)=e_{nu'}
    for nu in Ps:
        for mu in Ps:
            b=B[idx[nu],idx[mu]]
            if b!=0:
                de=sp.fraction(b)[1]; chk(sp.cancel(de.subs(s,1))!=0,f'b reg s=1 {nu},{mu}')
                chk(sp.cancel(b.subs(s,1))==(1 if mu==transpose(nu) else 0),f'P(1)=e {nu},{mu}')
    # Kostka-Foulkes K_{rho nu}(T): s_rho = sum_nu K P_nu(x;T);  s_rho = P(q=T,T), HL = P(q=0,T)
    Sm=sp.Matrix([[sp.cancel(sp.sympify(P[r].get(k,0)).subs(q,T)) for k in Ps] for r in Ps])
    Hm=sp.Matrix([[sp.cancel(sp.sympify(P[r].get(k,0)).subs(q,0)) for k in Ps] for r in Ps])
    KF=(Sm*Hm.inv()).applyfunc(sp.cancel)
    Tn=lambda l: t**nn(l)*s**nn(transpose(l))
    for lam in Ps:
        # (0) e* via (N)
        cN=[sp.cancel(t**(-nn(transpose(lam)))*sum(A[idx[lam],j]*Tn(Ps[j])*B[j,i] for j in range(len(Ps)))) for i in range(len(Ps))]
        if n<=5:
            v={():sp.Integer(1)}
            for k in reversed(lam): v=apply(k,v)
            cH={nu:sp.cancel(c*s**nn(nu)) for nu,c in v.items()}
            chk(all(sp.cancel(cN[idx[mu]]-cH.get(mu,0))==0 for mu in Ps),f'(N) vs Hikita {lam}')
        for mu in Ps:
            c=cN[idx[mu]]
            if dominates(mu,lam):
                chk(val(c,s)==nn(mu),f'val {lam},{mu}')
                d=lowcoef(c,s,nn(mu))
                aHL=sp.cancel(A[idx[lam],idx[transpose(mu)]].subs(s,0))
                chk(sp.cancel(d-t**(nn(transpose(mu))-nn(transpose(lam)))*aHL)==0,f'd=HL {lam},{mu}')
                # (4) tau=0 (t=oo) value of aHL = Kostka K_{mu,lam}
                aHLtau=sp.cancel(aHL.subs(t,1/T)); k0=sp.cancel(aHLtau.subs(T,0))
                D=nn(transpose(mu))-nn(transpose(lam))
                pred=sp.cancel(t**D*sum(KF[idx[r],idx[transpose(mu)]].subs(T,1/t)*KF[idx[transpose(r)],idx[lam]].subs(T,1) for r in Ps))
                chk(sp.cancel(d-pred)==0,f'd=Kostka formula {lam},{mu}')
                dp=sp.Poly(sp.cancel(d),t)
                chk(all(cc>=0 for cc in dp.coeffs()) and dp.eval(0)==1 and dp.eval(1)==em[lam].get(transpose(mu),0),f'd in N[t], d(0)=1, d(1)=M {lam},{mu}: {d}')
                chk(k0!=0 and k0>0,f'kostka {lam},{mu} {k0}')
                chk(sp.Poly(sp.cancel(d),t).degree()>=0 and sp.fraction(sp.cancel(d))[1].free_symbols==set() ,f'd poly {lam},{mu}: {d}')
            else:
                chk(c==0,f'supp {lam},{mu}')
            chk(sp.cancel(c.subs(s,1))==(1 if mu==lam else 0),f'V {lam},{mu}')
    print('n=',n,'OK',ok,'BAD',bad,flush=True)
print('TOTAL OK',ok,'BAD',bad)
