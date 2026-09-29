"""Day 209: coefficient form of (★2), n<=8, all i+j<=n, exact random rationals.
(1-t^n)V^n_ij = sum_p (s^2 t^{-i-j})^p [ V^{n-1-p}_ij (kg (t^i z)^{-p-1} + kd (t^j w)^{-p-1})
   - t^p V^{n-1-p}_{i-1,j} kg(i-1,j) (t^{i-1}z)^{-p-1} - t^p V^{n-1-p}_{i,j-1} kd(i,j-1) (t^{j-1}w)^{-p-1} ]"""
from fractions import Fraction as Fr
import math, random
random.seed(209)
def run(s,t,z,w,NM=8):
    poch=lambda x,N: math.prod((1-x*t**i for i in range(N)),start=Fr(1))
    al=lambda J: Fr(0) if J<0 else math.prod((s-t**i for i in range(1,J+1)),start=Fr(1))/poch(t,J)
    c=lambda n,j: Fr(0) if (j<0 or j>n) else poch(s,n-j)/poch(t,n-j)*(al(j)-s*t**(n-j)*al(j-1))
    Nn=lambda n,j: t**(-n*j)*c(n,j)
    K=lambda i,j: math.prod(((t**p*z-s*w)/(t**p*z-t**j*w) for p in range(i)),start=Fr(1))*math.prod(((s*z-t**r*w)/(t**i*z-t**r*w) for r in range(j)),start=Fr(1))
    def V(n,i,j):
        if i<0 or j<0 or n<0: return Fr(0)
        return K(i,j)*sum((s**(n1-i+n-n1-j)*t**(-(n1-i)*j-(n-n1-j)*i)*Nn(n1,i)*Nn(n-n1,j)*z**(-n1)*w**(-(n-n1)) for n1 in range(n+1)),Fr(0))
    kg=lambda i,j:(t**i*z-s*z)*(t**i*z-s*w)/(t**i*z*(t**i*z-t**j*w))
    kd=lambda i,j:(t**j*w-s*z)*(t**j*w-s*w)/(t**j*w*(t**j*w-t**i*z))
    ok=True
    for n in range(1,NM+1):
        for i in range(n+1):
            for j in range(n+1-i):
                r=Fr(0)
                for p in range(n):
                    a=(s**2*t**(-i-j))**p
                    r+=a*V(n-1-p,i,j)*(kg(i,j)*(t**i*z)**(-p-1)+kd(i,j)*(t**j*w)**(-p-1))
                    if i>0: r-=a*t**p*V(n-1-p,i-1,j)*kg(i-1,j)*(t**(i-1)*z)**(-p-1)
                    if j>0: r-=a*t**p*V(n-1-p,i,j-1)*kd(i,j-1)*(t**(j-1)*w)**(-p-1)
                ok&=(r==(1-t**n)*V(n,i,j))
    return ok
res=[run(Fr(random.randint(2,9),11),Fr(random.randint(2,9),13),Fr(random.randint(1,20),7),Fr(random.randint(21,40),9)) for _ in range(3)]
print(res,'ALL OK' if all(res) else 'FAIL')
