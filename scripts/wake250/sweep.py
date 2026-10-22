# Crown 1 kill test: F_k(M) =? C_k * w_M^a * sum_{nu ~ k} coef(k,nu) Psi_{nu,0}(M; base)
# base in {w, 1/w}; orientation up (nu>=k) / down (nu<=k); coef = Mobius mu or all +1 (control)
import sys, random
from fractions import Fraction as Fr
sys.path.insert(0,'../prove249c')
from eng import F, Psi, leq, subsets, cc
n=int(sys.argv[1]); spec=sys.argv[2]
if spec=='primes': w={i:Fr([5,7,11,13,17][i-1]) for i in range(1,n+1)}
else:
    random.seed(int(spec)); w={i:Fr(random.choice([5,7,11,13,17,19,23,29,31]),random.choice([37,41,43,47,53,59])) for i in range(1,n+1)}
wi={i:1/w[i] for i in w}
S=[tuple(s) for s in subsets(n)]
le=lambda a,b: leq(list(a),list(b))
def mobius(up):
    mu={}
    for a in S:
        for b in sorted([b for b in S if (le(a,b) if up else le(b,a))],key=lambda b:(sum(b),len(b)),reverse=not up):
            if a==b: mu[a,b]=1
            else:
                inn=[c for c in S if (a,c) in mu and c!=b and ((le(c,b)) if up else le(b,c))]
                mu[a,b]=-sum(mu[a,c] for c in inn)
    return mu
MU={True:mobius(True),False:mobius(False)}
def wp(M,a):
    p=Fr(1)
    for m in M: p*=w[m]**a
    return p
Ftab={(k,M):F(list(k),list(M),w) for k in S for M in S}
Ptab={(b,v,M):Psi(list(v),[0]*len(v),list(M),w if b=='w' else wi) for b in ('w','1/w') for v in S for M in S}
def test(base,up,a,plus,shift=None):
    """returns (#kappa with constant ratio, total)"""
    nc=0; tot=0
    for k in S:
        kk=k if shift is None else shift[k]
        rats=set(); bad=False
        for M in S:
            lhs=Ftab[k,M]
            rhs=wp(M,a)*sum((1 if plus else MU[up][kk,v])*Ptab[base,v,M] for v in S if (kk,v) in MU[up])
            if lhs==0 and rhs==0: continue
            if rhs==0 or lhs==0: bad=True; break
            rats.add(lhs/rhs)
        tot+=1
        if not bad and len(rats)==1:
            nc+=1
            pass  # cc-check done in ratios() below
    return nc,tot
def ratios(base,up,a):
    out={}
    for k in S:
        r=set()
        for M in S:
            lhs=Ftab[k,M]; rhs=wp(M,a)*sum(MU[up][k,v]*Ptab[base,v,M] for v in S if (k,v) in MU[up])
            if lhs or rhs: r.add(lhs/rhs if rhs else 'inf')
        out[k]=r
    return out
print('n=%d spec=%s'%(n,spec))
print('base orient twist(w_M^a) | Mobius: #const/total | all-plus control: #const/total')
win=[]
for base in ('w','1/w'):
    for up in (True,False):
        for a in (-2,-1,0,1,2):
            m=test(base,up,a,False); p=test(base,up,a,True)
            print('%-4s %-5s a=%+d | %d/%d | %d/%d'%(base,'up' if up else 'down',a,m[0],m[1],p[0],p[1]))
            if m[0]==m[1]: win.append((base,up,a))
print('winning combos:',win)
for base,up,a in win:
    R=ratios(base,up,a)
    ok=all(len(R[k])==1 and list(R[k])[0]==cc(list(k),w) for k in S)
    print('combo',base,up,a,': ratio == cc(kappa)=(-1)^|k| prod_{k in kappa} w_1..w_k for all kappa:',ok)
    # perturbation control: compare F_k with xi^{k'} for k' != k (any permutation-free: all pairs)
    npass=0; ntot=0
    for k in S:
        for k2 in S:
            if k2==k: continue
            ntot+=1; rats=set(); bad=False
            for M in S:
                lhs=Ftab[k,M]; rhs=wp(M,a)*sum(MU[up][k2,v]*Ptab[base,v,M] for v in S if (k2,v) in MU[up])
                if lhs==0 and rhs==0: continue
                if lhs==0 or rhs==0: bad=True; break
                rats.add(lhs/rhs)
            if not bad and len(rats)==1: npass+=1
    print('  perturbation control (F_k vs xi^{k\'}, k\'!=k): spuriously constant %d/%d'%(npass,ntot))
