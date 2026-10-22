# Independent K_T Schubert classes of LG(n,2n)=Sp_2n/P (P omits alpha_n) via Graham-Willems localization.
# Roots of C_n as integer vectors; characters evaluated at t_p = e^{eps_p} (Fractions).
import itertools
from fractions import Fraction as Fr
def mul(a,b):  # signed perms as tuples of images of 1..n ; (a*b)(i)=a(b(i))
    return tuple((1 if x>0 else -1)*a[abs(x)-1] for x in b)
def simple(n,i):
    if i<n:
        l=list(range(1,n+1)); l[i-1],l[i]=l[i],l[i-1]; return tuple(l)
    l=list(range(1,n+1)); l[n-1]=-n; return tuple(l)
def act(w,vec):  # w acting on vector in eps basis
    n=len(w); out=[0]*n
    for i in range(n):
        out[abs(w[i])-1]+= (1 if w[i]>0 else -1)*vec[i]
    return tuple(out)
def posroots(n):
    R=[]
    for i in range(n):
        for j in range(i+1,n):
            v=[0]*n; v[i]=1; v[j]=-1; R.append(tuple(v))
            v=[0]*n; v[i]=1; v[j]=1; R.append(tuple(v))
        v=[0]*n; v[i]=2; R.append(tuple(v))
    return R
def ispos(v):
    for x in v:
        if x: return x>0
    return False
def alpha(n,i):
    v=[0]*n
    if i<n: v[i-1]=1; v[i]=-1
    else: v[n-1]=2
    return tuple(v)
class Cn:
    def __init__(s,n):
        s.n=n; s.R=posroots(n); s.s={i:simple(n,i) for i in range(1,n+1)}
        s.W=[tuple(sg*p for sg,p in zip(sgn,perm)) for perm in itertools.permutations(range(1,n+1)) for sgn in itertools.product([1,-1],repeat=n)]
        s.len={w:sum(1 for b in s.R if not ispos(act(w,b))) for w in s.W}
        # cosets: key = sign vector w(omega_n)
        s.rep={}
        for w in s.W:
            k=act(w,(1,)*n)
            if k not in s.rep or s.len[w]<s.len[s.rep[k]]: s.rep[k]=w
    def word(s,w):
        wd=[]
        while s.len[w]>0:
            for i in range(1,s.n+1):
                if not ispos(act(w,alpha(s.n,i))):
                    wd.append(i); w=mul(w,s.s[i]); break
        return wd[::-1]
    def hecke(s,x,i):
        y=mul(x,s.s[i]); return y if s.len[y]>s.len[x] else x
    def ech(s,v,t):  # e^v at t
        p=Fr(1)
        for vi,ti in zip(v,t): p*=ti**vi
        return p
    def willems(s,v,u,t):
        """[O_{X^v}]|_u , X^v opposite Schubert (codim l(v)), Graham-Willems: sum over subwords with Hecke product v"""
        wd=s.word(u); n=s.n
        betas=[]; pre=tuple(range(1,n+1))
        for i in wd:
            betas.append(act(pre,alpha(n,i))); pre=mul(pre,s.s[i])
        e=tuple(range(1,n+1)); tot=Fr(0)
        # DP over word: state = (hecke product, count) -> weight
        st={(e,0):Fr(1)}
        for j,i in enumerate(wd):
            new=dict(st)
            f=1-s.ech(betas[j],t)
            for (x,k),val in st.items():
                y=s.hecke(x,i); key=(y,k+1)
                new[key]=new.get(key,0)+val*f
            st=new
        lv=s.len[v]
        for (x,k),val in st.items():
            if x==v: tot+=(-1)**(k-lv)*val
        return tot
    def tangent(s,u):  # tangent weights at uP: u(-beta), beta in R+ \ R+_P
        RP=[b for b in s.R if b[0:s.n].count(0)==s.n-2 and sum(b)==0 and len([x for x in b if x])==2 and 1 in b and -1 in b]
        return [act(u,tuple(-x for x in b)) for b in s.R if b not in RP]
