import sys, itertools
sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/prove249d')
from dem import Di, word
sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/wake251/orderpoly')
from orderpoly import shifted, boxless, longest_chain
def chiw(kap,n,f):
    for i in reversed(word(set(kap),n)): f=Di(f,i,n)
    return sum(f.values())
n=int(sys.argv[1])
for r in range(1,n+1):
  for kap in itertools.combinations(range(1,n+1),r):
    if max(kap)!=n: continue
    lam=sorted(kap,reverse=True); rr=longest_chain(shifted(lam),boxless); rest=lam[1:]; rp=longest_chain(shifted(rest),boxless) if rest else 0
    Z=[]
    for k in range(n+1):
        zs=[a for a in range(1,n+2) if chiw(kap,n,{tuple(-a-s for s in S):1 for S in itertools.product((0,1),repeat=n) if sum(S)==k})==0]
        # first a>=1 not zero
        first=next(a for a in range(1,3*n) if chiw(kap,n,{tuple(-a-s for s in S):1 for S in itertools.product((0,1),repeat=n) if sum(S)==k})!=0)
        Z.append(first-1)
    print(kap,'r',rr,"r'",rp,'window(r)=',(rr+1)//2,"ceil(r'/2)=",(rp+1)//2,' maxzero a per k=0..n:',Z)
