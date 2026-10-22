# Task 2: full chain at fixed points, all classes from Graham-Willems (independent of Psi), AB localisation.
import sys, random
from fractions import Fraction as Fr
sys.path.insert(0,'../prove249c')
from eng import Psi, F, ell, cc, subsets, leq
from kschub import Cn, act
n=int(sys.argv[1]); seed=int(sys.argv[2]) if len(sys.argv)>2 else 0
if seed==0: w={i:Fr([5,7,11,13,17][i-1]) for i in range(1,n+1)}
else:
    random.seed(seed); w={i:Fr(random.choice([5,7,11,13,17,19,23,29,31]),random.choice([37,41,43,47,53,59])) for i in range(1,n+1)}
C=Cn(n); S=[tuple(x) for x in subsets(n)]
def sv(M):
    v=[1]*n
    for k in M: v[n-k]=-1
    return tuple(v)
comp=lambda M: tuple(k for k in range(1,n+1) if k not in M)
t=[w[n-p] for p in range(n)]          # e^{eps_p}=w_{n+1-p}  (the dictionary found in dfree1)
ti=[1/x for x in t]
Oup={(v,M):C.willems(C.rep[sv(v)],C.rep[sv(M)],t) for v in S for M in S}        # [O_{X^v}]|_M
# Schubert X_u = w0 X^{w0 u}; w0=-1: coset of sign vector -> complement; characters t->1/t
Odn={(u,M):C.willems(C.rep[sv(comp(u))],C.rep[sv(comp(M))],ti) for u in S for M in S}  # [O_{X_u}]|_M
def ech(v,tt):
    p=Fr(1)
    for a,b in zip(v,tt): p*=b**a
    return p
out=[]
def pr(*a):
    s=' '.join(str(x) for x in a); print(s); out.append(s)
pr('n=%d seed=%d'%(n,seed))
pr('X_u support = {M<=u}:',all((Odn[u,M]!=0)==leq(list(M),list(u)) for u in S for M in S))
best=None
for sgn in (1,-1):
    eul={M:Fr(1) for M in S}
    for M in S:
        for tau in C.tangent(C.rep[sv(M)]): eul[M]*=1-ech(tuple(sgn*x for x in tau),t)
    chi=lambda f: sum(f[M]/eul[M] for M in S)
    one=chi({M:1 for M in S})
    rich=all(chi({M:Odn[u,M]*Oup[v,M] for M in S})==(1 if leq(list(v),list(u)) else 0) for u in S for v in S)
    pr('AB sign %+d: chi(O)=%s  chi(O_{X_u}O_{X^v})=[v<=u] all u,v: %s'%(sgn,one,rich))
    if one==1 and rich: best=(sgn,eul)
sgn,eul=best
chi=lambda f: sum(f[M]/eul[M] for M in S)
# dual basis xi^v to {O_{X_u}}: solve sum_M O_u(M) xi^v(M)/eul(M) = delta_uv  (matrix inverse)
import sympy
A=sympy.Matrix([[sympy.Rational(Odn[u,M].numerator,Odn[u,M].denominator)/sympy.Rational(eul[M].numerator,eul[M].denominator) for M in S] for u in S])
X=A.inv()   # columns: xi^v as vector over M  (A*X=I => sum_M A[u,M] X[M,v]=delta)
xi={(v,M):Fr(int(sympy.fraction(X[i,j])[0]),int(sympy.fraction(X[i,j])[1])) for j,v in enumerate(S) for i,M in enumerate(S)}
# Mobius of dominance on S
le=lambda a,b: leq(list(a),list(b))
mu={}
for a in S:
    for b in sorted([b for b in S if le(a,b)],key=lambda b:(sum(b),len(b))):
        mu[a,b]=1 if a==b else -sum(mu[a,c] for c in S if (a,c) in mu and c!=b and le(c,b))
mob_ok=all(xi[v,M]==sum(mu[v,x]*Oup[x,M] for x in S if (v,x) in mu) for v in S for M in S)
pr('xi^v (chi_T-dual of O_{X_u}) == sum_{x>=v} mu(v,x) O_{X^x}  :',mob_ok)
pr('mu values on W^P intervals:',sorted(set(mu.values())), ' mu==(-1)^{len diff} where nonzero:',all(m==0 or m==(-1)**(sum(b)-sum(a)) for (a,b),m in mu.items()))
# Psi dictionary (independent check repeated)
pr('Psi_{nu,0}(M;1/w)==O_{X^nu}|_M all:',all(Psi(list(v),[0]*len(v),list(M),{k:1/w[k] for k in w})==Oup[v,M] for v in S for M in S))
pr('Psi_{nu,0}(M;w) ==O_{X^nu}|_M all (neg control):',all(Psi(list(v),[0]*len(v),list(M),w)==Oup[v,M] for v in S for M in S),
   ' #mismatch=',sum(Psi(list(v),[0]*len(v),list(M),w)!=Oup[v,M] for v in S for M in S),'/',len(S)**2)
# Missing lemma: w_M^2 F_k(M) = cc(k) xi^k(M)
wM=lambda M: eval('*'.join(['1']+['w[%d]'%m for m in M]))
lem=[(k,M) for k in S for M in S if wM(M)**2*F(list(k),list(M),w)!=cc(list(k),w)*xi[k,M]]
pr('LEMMA  w_M^2 F_k(M) == cc(k) xi^k(M)  mismatches: %d/%d'%(len(lem),len(S)**2))
# full chain: chi_T(O_{X_k} * w^2 F_l)/cc(l) = delta  (Theorem-D form of (D))
ch=[(k,l) for k in S for l in S if chi({M:Odn[k,M]*wM(M)**2*F(list(l),list(M),w) for M in S})!=(cc(list(l),w) if k==l else 0)]
pr('chi_T(O_{X_k} . w^2 F_l) == cc(k) delta_{kl}  mismatches: %d/%d'%(len(ch),len(S)**2))
# Motegi explicit (D) functional vs geometric functional: ell(k,L) =? O_{X_k}(L) w_L^2 /(cc(k) eul(L))
ed=[(k,L) for k in S for L in S if ell(list(k),list(L),n,w)!=Odn[k,L]*wM(L)**2/(cc(list(k),w)*eul[L])]
pr('explicit ell_k(e_L) == O_{X_k}|_L w_L^2/(cc(k) eul_L)  mismatches: %d/%d'%(len(ed),len(S)**2))
# split of the explicit formula: (i) Schubert-class part: Psi_{k^c,0}(L^c) vs O_{X_k}(L)
w_inv={k:1/w[k] for k in w}
r1=set()
for k in S:
    for L in S:
        a=Psi(list(comp(k)),[0]*len(comp(k)),list(comp(L)),w); b=Odn[k,L]
        if a or b: r1.add((k, a/b if b else 'inf'))
from collections import defaultdict
d=defaultdict(set)
for k,r in r1: d[k].add(r)
pr('(i) Psi_{k^c,0}(L^c; w) / O_{X_k}|_L constant in L for every k:',all(len(v)==1 for v in d.values()),' constants==1:',all(v=={1} for v in d.values()))
d2=defaultdict(set)
for k in S:
    for L in S:
        a=Psi(list(comp(k)),[0]*len(comp(k)),list(comp(L)),w_inv); b=Odn[k,L]
        if a or b: d2[k].add(a/b if b else 'inf')
pr('(i-ctrl) same with base 1/w constant per k:',sum(len(v)==1 for v in d2.values()),'/',len(S))
# (ii) L-only identity: cc(L)/(F_L(L) Psi_{L^c,0}(L^c)) == w_L^2/eul(L)
ii=[L for L in S if cc(list(L),w)/(F(list(L),list(L),w)*Psi(list(comp(L)),[0]*len(comp(L)),list(comp(L)),w))!=wM(L)**2/eul[L]]
pr('(ii) cc(L)/(F_L(L)Psi_{L^c,0}(L^c)) == w_L^2/eul_L  mismatches: %d/%d'%(len(ii),len(S)))
open('dfree_out_n%d_s%d.txt'%(n,seed),'w').write('\n'.join(out)+'\n')
