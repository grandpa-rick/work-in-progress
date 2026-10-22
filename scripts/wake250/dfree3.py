# Task 2b: do the independently computed ideal-sheaf classes cc(l) xi^l satisfy Lemmas A_m / B (left W-action) -- the
# recursion proved for g_l = w^2 F_l?  xi from Graham-Willems + AB duality only (no Psi, no F).
import sys
from fractions import Fraction as Fr
sys.path.insert(0,'../prove249c')
from eng import cc, subsets, F
from kschub import Cn, act
import sympy
n=int(sys.argv[1]); C=Cn(n); S=[tuple(x) for x in subsets(n)]
def sv(M):
    v=[1]*n
    for k in M: v[n-k]=-1
    return tuple(v)
comp=lambda M: tuple(k for k in range(1,n+1) if k not in M)
def ech(v,tt):
    p=Fr(1)
    for a,b in zip(v,tt): p*=b**a
    return p
def xi_at(w):
    t=[w[n-p] for p in range(n)]; ti=[1/x for x in t]
    Odn={(u,M):C.willems(C.rep[sv(comp(u))],C.rep[sv(comp(M))],ti) for u in S for M in S}
    eul={}
    for M in S:
        e=Fr(1)
        for tau in C.tangent(C.rep[sv(M)]): e*=1-ech(tau,t)
        eul[M]=e
    A=sympy.Matrix([[sympy.Rational(Odn[u,M].numerator,Odn[u,M].denominator)/sympy.Rational(eul[M].numerator,eul[M].denominator) for M in S] for u in S])
    X=A.inv()
    return {(v,M):cc(list(v),w)*Fr(int(sympy.fraction(X[i,j])[0]),int(sympy.fraction(X[i,j])[1])) for j,v in enumerate(S) for i,M in enumerate(S)}
w={i:Fr([5,7,11,13,17][i-1]) for i in range(1,n+1)}
G=xi_at(w)   # G[l,M] = cc(l) xi^l(M)   (claimed = g_l = w^2 F_l)
out=[]
def pr(*a):
    s=' '.join(str(x) for x in a); print(s); out.append(s)
pr('n=%d  G_l := cc(l) xi^l (geometric, Psi/F-free)'%n)
fails=0; tot=0
for m in range(1,n):
    ws=dict(w); ws[m],ws[m+1]=w[m+1],w[m]; Gs=xi_at(ws)
    sw=lambda L: tuple(sorted((m+1 if x==m else m if x==m+1 else x) for x in L))
    for l in S:
        lhs={L:Gs[l,sw(L)] for L in S}   # (s.G_l)(L)=s(G_l(sL))
        if m+1 not in l: rhs={L:G[l,L] for L in S}; case='a'
        elif m not in l:
            lp=tuple(sorted(m if x==m+1 else x for x in l)); rhs={L:G[l,L]+(w[m]-w[m+1])*G[lp,L] for L in S}; case='b'
        else: rhs={L:w[m+1]/w[m]*G[l,L] for L in S}; case='c'
        ok=lhs==rhs; tot+=1; fails+=not ok
        if not ok: pr('  FAIL A_%d(%s) l=%s'%(m,case,l))
ws=dict(w); ws[1]=1/w[1]; Gs=xi_at(ws); z=w[1]
sg=lambda L: tuple(sorted(set(L)^{1}))
for l in S:
    r=len(l); lhs={L:Gs[l,sg(L)] for L in S}
    if 1 not in l: rhs={L:z**(-2*(r+1))*G[l,L] for L in S}; case='a'
    else:
        lm=tuple(x for x in l if x!=1); rhs={L:z**(-2*r)*G[l,L]-(z*z-1)*z**(-2*r-1)*G[lm,L] for L in S}; case='b'
    ok=lhs==rhs; tot+=1; fails+=not ok
    if not ok: pr('  FAIL B(%s) l=%s'%(case,l))
pr('Lemmas A_m(a,b,c)+B(a,b) hold for geometric cc*xi: %d/%d instances pass'%(tot-fails,tot))
# negative control: plain xi (without cc normalisation) should fail some
Gn={k:v/cc(list(k[0]),w) for k,v in G.items()}
pr('top class: G_[n] supported only at M=[n]:',all(G[tuple(range(1,n+1)),M]==0 for M in S if M!=tuple(range(1,n+1))))
wM=lambda M: eval('*'.join(['1']+['w[%d]'%x for x in M]))
pr('G_[n]([n]) == w^2 F_[n]([n]):',G[tuple(range(1,n+1)),tuple(range(1,n+1))]==wM(range(1,n+1))**2*F(list(range(1,n+1)),list(range(1,n+1)),w))
open('dfree3_n%d.txt'%n,'w').write('\n'.join(out)+'\n')
