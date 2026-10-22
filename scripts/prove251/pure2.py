import sys, itertools
sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/prove249d')
from dem import Di
from span import G, Fb
for a,ks in [(1,[1,2,3]),(2,[3,4,5]),(3,[5,6,7]),(4,[7,8])]:
    for k in ks:
        g=G(k,a); pure={tuple([-a]*k):1}
        target = pure if k>=2*a else {**{t:-v for t,v in Fb(k,a-1).items()}, tuple([-a]*k):1+(-Fb(k,a-1).get(tuple([-a]*k),0))*0}
        if k==2*a-1:
            t=dict(pure)
            for kk,v in Fb(k,a-1).items(): t[kk]=t.get(kk,0)-v
            target={kk:v for kk,v in t.items() if v}
        print('a',a,'k',k,'claim', 'pure' if k>=2*a else 'e^{-aw}-F_{a-1}', 'OK' if g==target else 'FAIL')
