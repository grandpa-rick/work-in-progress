import sys, random
from fractions import Fraction as Fr
sys.path.insert(0,'../prove249c')
from eng import Cfun, subsets
n=3; random.seed(5)
w={i:Fr(random.randint(2,30),random.randint(2,30)) for i in range(1,n+1)}
mw={i:-w[i] for i in w}
for lam,mu,kap in [([],[],[3]),([],[],[1,2]),([1],[],[2,3]),([1],[2],[1,3]),([],[],[1,2,3])]:
    print(lam,mu,kap,[ (Cfun(lam,mu,kap,c,n,mw)/Cfun(lam,mu,kap,c,n,w)) if Cfun(lam,mu,kap,c,n,w)!=0 else 'z' for c in range(-3,4)])
