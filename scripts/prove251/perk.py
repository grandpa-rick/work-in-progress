import sys, itertools
sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/prove249d')
from dem import Di, word
sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/wake251/orderpoly')
from orderpoly import shifted, boxless, longest_chain
def chiw(kap,n,f):
    ws=word(set(kap),n)
    for i in reversed(ws): f=Di(f,i,n)
    return sum(f.values())
n=int(sys.argv[1])
for r in range(1,n+1):
  for kap in itertools.combinations(range(1,n+1),r):
    rest=sorted(kap,reverse=True)[1:]; rp=longest_chain(shifted(rest),boxless) if rest else 0
    rows=[]
    for a in range(0,n+2):
        f={}
        for S in itertools.product((0,1),repeat=n): f[tuple(-a-s for s in S)]=f.get(tuple(-a-s for s in S),0)+1
        perk=[]
        for k in range(n+1):
            g={tuple(-a-s for s in S):1 for S in itertools.product((0,1),repeat=n) if sum(S)==k}
            perk.append(chiw(kap,n,g))
        rows.append(perk)
    print(kap,"r'=",rp, ' a=0..: chi(wedge^k U(-a)) k=0..n:', rows)
