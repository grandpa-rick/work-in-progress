from psi import *
import sys
def lead0(expr,u):
    e=sp.cancel(sp.together(expr)); num,den=sp.fraction(e)
    pn=sp.Poly(num,u); pd=sp.Poly(den,u)
    on=min(m[0] for m in pn.monoms()); od=min(m[0] for m in pd.monoms())
    cn=sum(c*u**0 for (m,),c in zip(pn.monoms(),pn.coeffs()) if m==on)
    cd=sum(c for (m,),c in zip(pd.monoms(),pd.coeffs()) if m==od)
    return on-od, sp.factor(cn/cd)
def B(s,L,M,th=None):
    # B_L via Q formula: sum_j a1..a_{j-1} prod_{i>=j}^{L}(x_i-a_i), labels (s..s+L-1)
    kap=tuple(range(s,s+L)); a=[w(s+i) for i in range(L)]
    tot=0
    for h in itertools.permutations(M,L):
        if any(h[j]<kap[j] for j in range(L)): continue
        x=[w(m) for m in h]
        Q=0
        for j in range(L+1):
            t=sp.Mul(*a[:j])
            for i in range(j,L): t*=(x[i]-a[i])
            Q+=t
        tot+=T0(kap,h,M)*Q
    return sp.factor(sp.cancel(sp.together(tot)))
if __name__=="__main__":
  s=int(sys.argv[1]); L=int(sys.argv[2]); M=tuple(map(int,sys.argv[3].split(','))); n=int(sys.argv[4])
  o,c=lead0(Psi(tuple(range(s,s+L+1)),(0,)*(L+1),M+(n,)),w(n))
  b=B(s,L,M)
  print('order',o); print('lead',c); print('B',b); print('ratio',sp.factor(c/b))
