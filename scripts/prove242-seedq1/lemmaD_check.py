# Day 242 cycle 2: Lemma D without Lusztig 1.2.3. Form B(x,F_w):=kappa_w c_w(x), c_w from ITERATED r (braided tensor algebra hom),
# computed independently of the shuffle code. Checks: (1) c_w(F_u) = shuffle coeff (Lemma A dictionary), (2) B symmetric,
# (3) B(x, y'y'') = (r x, y'⊗y'') factorwise, (4) rank Gram = Kostant dim.
import sys; sys.path.insert(0,'code/seedq1')
from fractions import Fraction as Fr
from itertools import permutations
from collections import defaultdict
import sympy as sp
exec(open('code/seedq1/shuffle_checks.py').read().split('types = {')[0])  # shmul/lin/theta_word only
def deg(w, n): 
    v=[0]*n
    for x in w: v[x]+=1
    return v
def pair(a,b,C): return sum(a[i]*C[i][j]*b[j] for i in range(len(a)) for j in range(len(a)))
def r_word(w, C, qq, n):  # r(F_w) in 'f⊗'f, braided: (a⊗b)(c⊗d)=qq^{(|b|,|c|)} ac⊗bd
    X={((),()):1}
    for x in w:
        Y=defaultdict(lambda:0)
        for (a,b),c in X.items():
            Y[(a+(x,),b)] += c*qq**pair(deg(b,n),deg((x,),n),C)   # (a⊗b)(F_x⊗1)
            Y[(a,b+(x,))] += c                                     # (a⊗b)(1⊗F_x)
        X=Y
    return X
def cw_free(u, w, C, qq, n):  # c_w(F_u) via iterated r: peel first letter of w with r, recurse
    if len(w)==0: return 1 if len(u)==0 else 0
    if len(w)==1: return 1 if u==w else 0
    tot=0
    for (a,b),c in r_word(u,C,qq,n).items():
        if a==(w[0],): tot += c*cw_free(b,w[1:],C,qq,n)
    return tot
types={"A3":([[2,-1,0],[-1,2,-1],[0,-1,2]],[1,1,1]),"B2":([[2,-2],[-1,2]],[1,2]),"G2":([[2,-3],[-1,2]],[1,3])}
qq=Fr(3,7); bad=0; tests=0
for name,(A,d) in types.items():
    n=len(A); C=[[d[i]*A[i][j] for j in range(n)] for i in range(n)]
    kap=[1/(1-qq**(-2*d[i])) for i in range(n)]
    contents={"A3":[(1,2,1),(2,2,1),(1,2,2),(2,1,1)],"B2":[(2,1),(1,2),(3,1),(2,2)],"G2":[(2,1),(3,1),(4,1),(3,2)]}[name]
    for ct in contents:
        base=[i for i in range(n) for _ in range(ct[i])]
        W=sorted(set(permutations(base)))
        kw=1
        for i in range(n): kw*=kap[i]**ct[i]
        G=[[kw*cw_free(u,w,C,qq,n) for w in W] for u in W]
        # (1) dictionary with shuffle
        for u in W:
            th=theta_word(u,C,qq)
            for w in W:
                tests+=1
                if th.get(w,0)!=cw_free(u,w,C,qq,n): bad+=1
        sym=all(G[i][j]==G[j][i] for i in range(len(W)) for j in range(len(W)))
        # (3) Hopf property: B(F_u, F_a F_b) = sum r(F_u)_{(x1,x2)} B(x1,F_a)B(x2,F_b)
        hop=True
        for u in W:
            ru=r_word(u,C,qq,n)
            for w in W:
                for m in range(1,len(w)):
                    a,b=w[:m],w[m:]
                    ka=1;kb=1
                    for x in a: ka*=kap[x]
                    for x in b: kb*=kap[x]
                    rhs=sum(c*ka*cw_free(x1,a,C,qq,n)*kb*cw_free(x2,b,C,qq,n) for (x1,x2),c in ru.items() if len(x1)==m)
                    if rhs!=kw*cw_free(u,w,C,qq,n): hop=False
        rk=sp.Matrix(G).rank()
        print(name,ct,"#words",len(W),"symmetric",sym,"hopf",hop,"rank",rk)
print("dictionary c_w(r) == shuffle coeff:",tests-bad,"/",tests)
# CONTROL: non-symmetric braiding matrix (A2 Cartan a_ij itself made asymmetric: B2 unsymmetrised) must break symmetry
Cn=[[2,-2],[-1,2]]; W=sorted(set(permutations((0,0,1))))
Gn=[[cw_free(u,w,Cn,qq,2) for w in W] for u in W]
print("CONTROL unsymmetrised B2 Cartan as braiding: symmetric?", all(Gn[i][j]==Gn[j][i] for i in range(3) for j in range(3)))
