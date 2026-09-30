"""Day 212: concrete-I checks at exact rational points (Fractions), all I with small |I|:
(P5) V_I^{(n)} = [x^n] K_I s^{-(l-1)S} t^{S^2-Σi^2} ∏ C_{i_c}(X_c), X_c = s^{l-1}t^{-S}x/z_c
(P6) κ'_c K_{I-e_c}/K_I = (A_c-st)/A_c ∏_{c'≠c} A_{c'}(A_cρ-1)/(A_cρ-A_{c'}), ρ=z_c/z_{c'}
(★ℓ) coefficient form, n<=NN (end-to-end for §3-§4, independent of normalisation algebra)."""
import itertools
from fractions import Fraction as Fr
def helpers(ell,s,t,Z):
    def poch(x,N):
        r=Fr(1)
        for i in range(N): r*=1-x*t**i
        return r
    def al(J):
        if J<0: return Fr(0)
        r=Fr(1)
        for i in range(1,J+1): r*=s-t**i
        return r/poch(t,J)
    def cf(n,j): return Fr(0) if (j<0 or j>n) else poch(s,n-j)/poch(t,n-j)*(al(j)-s*t**(n-j)*al(j-1))
    def N(n,j): return t**(-n*j)*cf(n,j)
    def Kp(i,j,z,w):
        r=Fr(1)
        for p in range(i): r*=(t**p*z-s*w)/(t**p*z-t**j*w)
        for q in range(j): r*=(s*z-t**q*w)/(t**i*z-t**q*w)
        return r
    def KI(I):
        r=Fr(1)
        for a in range(ell):
            for b in range(a+1,ell): r*=Kp(I[a],I[b],Z[a],Z[b])
        return r
    def comps(n,l):
        if l==1: yield (n,); return
        for a in range(n+1):
            for r in comps(n-a,l-1): yield (a,)+r
    def V(n,I):
        if min(I)<0: return Fr(0)
        tot=Fr(0)
        for ns in comps(n,ell):
            if any(ns[c]<I[c] for c in range(ell)): continue
            v=s**((ell-1)*sum(ns[c]-I[c] for c in range(ell)))*t**(-sum((ns[c]-I[c])*I[d] for c in range(ell) for d in range(ell) if d!=c))
            for c in range(ell): v*=N(ns[c],I[c])*Z[c]**(-ns[c])
            tot+=v
        return KI(I)*tot
    def kap(I,c):
        gam=[t**I[d]*Z[d] for d in range(ell)]
        num=Fr(1); den=gam[c]
        for d in range(ell): num*=gam[c]-s*Z[d]
        for d in range(ell):
            if d!=c: den*=gam[c]-gam[d]
        return num/den
    return cf,KI,V,kap
ok=True
def rep(msg,v):
    global ok; ok&=v; print(msg,'OK' if v else 'FAIL',flush=True)
pts=[(Fr(3,7),Fr(5,11),[Fr(2,3),Fr(9,8),Fr(13,5),Fr(17,19)]),(Fr(-2,9),Fr(7,4),[Fr(5,2),Fr(-3,7),Fr(11,6),Fr(1,13)])]
NN=6
for pi,(s,t,Zall) in enumerate(pts):
    for ell in (1,2,3,4):
        Z=Zall[:ell]; cf,KI,V,kap=helpers(ell,s,t,Z)
        Is=[I for I in itertools.product(range(4),repeat=ell) if sum(I)<=(4 if ell<4 else 3)]
        p5=p6=st=True
        for I in Is:
            S=sum(I); A=[t**i for i in I]; cc=s**ell*t**(-S)
            r=[cc/(s*Z[c]) for c in range(ell)]
            poly=[Fr(1)]+[Fr(0)]*NN
            for c in range(ell):
                ser=[cf(n,I[c])*r[c]**n for n in range(NN+1)]
                poly=[sum(poly[a]*ser[n-a] for a in range(n+1)) for n in range(NN+1)]
            pref=KI(I)*s**(-(ell-1)*S)*t**(S*S-sum(i*i for i in I))
            for n in range(NN+1):
                if pref*poly[n]!=V(n,I): p5=False
            for c in range(ell):
                if I[c]==0: continue
                J=list(I); J[c]-=1; J=tuple(J)
                lhs=kap(J,c)*KI(J)/KI(I)
                rhs=(A[c]-s*t)/A[c]
                for d in range(ell):
                    if d!=c:
                        rho=Z[c]/Z[d]; rhs*=A[d]*(A[c]*rho-1)/(A[c]*rho-A[d])
                if lhs!=rhs: p6=False
            gam=[A[c]*Z[c] for c in range(ell)]; ks=[kap(I,c) for c in range(ell)]
            for n in range(NN+1):
                tot=Fr(0)
                for p in range(n):
                    tot+=cc**p*V(n-1-p,I)*sum(ks[c]*gam[c]**(-p-1) for c in range(ell))
                    for c in range(ell):
                        if I[c]==0: continue
                        J=list(I); J[c]-=1; J=tuple(J)
                        tot-=cc**p*t**p*V(n-1-p,J)*kap(J,c)*(t**(I[c]-1)*Z[c])**(-p-1)
                if tot!=(1-t**n)*V(n,I): st=False
        rep(f'pt{pi} ell={ell} ({len(Is)} I, n<={NN}) (P5)',p5)
        rep(f'pt{pi} ell={ell} (P6)',p6)
        rep(f'pt{pi} ell={ell} (★ℓ coeff form)',st)
print('ALL OK' if ok else 'FAIL')
