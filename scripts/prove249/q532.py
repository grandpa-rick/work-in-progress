# Q532 checks: Clio's slope -1 collapse, breaking controls, u->oo peel of consecutive F-blocks
import sympy as sp, itertools, random
from fractions import Fraction as Fr
random.seed(532)
def T0v(kap,h,M,x):
    r=len(kap); X=[x[m] for m in h]; e=1
    for j in range(r):
        e*=(1-X[j]**-2)
        for k in range(1,kap[j]): e*=(1-x[k]/X[j])
    for j in range(r):
        for k in range(j+1,r): e*=(1-1/(X[j]*X[k]))/(1-X[j]/X[k])
    for t in M:
        if t in h: continue
        for j in range(r): e*=(1-1/(x[t]*X[j]))/(1-x[t]/X[j])
    return e
def Psi(kap,th,M,x):
    s=0
    for h in itertools.permutations(M,len(kap)):
        if any(h[j]<kap[j] for j in range(len(kap))): continue
        X=[x[m] for m in h]; term=T0v(kap,h,M,x)
        for j in range(len(kap)): term*=X[j]**th[j]
        s+=term
    return s
def pt(N): return {k:Fr(random.randint(2,97),random.randint(2,97))+k for k in range(1,N+1)}
print("== (1) identity Psi_{kappa, slope -1 on block} = C * Psi_{kappa, const e on block} ==")
cases=[ # kappa, block positions (0-based p..q), t0, other twists, M
 ((2,3,4),(0,2),1,{}, (2,3,4,5,6)),
 ((3,4),(0,1),2,{}, (1,3,4,5)),
 ((1,3,4,5),(1,3),-1,{0:2},(1,2,3,4,5,6)),
 ((2,3,4,5),(0,3),0,{},(2,3,4,5,6)),
 ((1,2,5,6),(2,3),3,{0:-1,1:1},(1,2,4,5,6)),
]
for kap,(p,q),t0,oth,M in cases:
    for trial in range(2):
        x=pt(max(M)+1)
        th1=[oth.get(j,0) for j in range(len(kap))]; th2=list(th1)
        b=q-p+1; e=t0-(b-1)
        for j in range(p,q+1): th1[j]=t0-(j-p); th2[j]=e
        a=[x[kap[p]+m] for m in range(b-1)]
        C=1
        for m in range(b-1): C*=a[m]**(b-1-m)
        L=Psi(kap,th1,M,x); R=C*Psi(kap,th2,M,x)
        print(kap,'block',(p,q),'t0',t0,'M',M,'diff==0:',L-R==0)
print("== breaking controls (must be nonzero) ==")
for kap,th1,th2,M,lab in [((2,3,4),(1,0,-1),(-1,-1,-1),(2,3,4,5,6),'C omitted (C:=1)'),
                          ((2,3,4),(2,0,-2),(-2,-2,-2),(2,3,4,5,6),'slope -2'),
                          ((2,4,5),(1,0,-1),(-1,-1,-1),(2,3,4,5,6),'non-consecutive labels')]:
    x=pt(7)
    if lab=='C omitted (C:=1)': C=1
    else: C=x[kap[0]]**2*x[kap[0]+1]
    print(lab, 'diff!=0:', Psi(kap,th1,M,x)-C*Psi(kap,th2,M,x)!=0)
