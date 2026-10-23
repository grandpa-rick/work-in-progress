import sys; sys.path.insert(0,'.')
N=int(sys.argv[1]); sys.argv=[sys.argv[0],str(N)]
import itertools, collections, gate6 as G, brpieri as B
from fractions import Fraction as Fr
st=collections.Counter(); ex=[]
for L in range(2,N+1):
  for c in itertools.combinations(range(1,N+1),L):
    k=tuple(sorted(c,reverse=True)); rp=G.rp(k); Ms=2*((rp+1)//2)+1; K=B.pieri_dim(k,N)
    cs=all(v*(-1)**(sum(k)-sum(r)-1)>0 for r,v in K.items())
    terms={r:v*G.H(r,Fr(-Ms,2)) for r,v in K.items()}
    nz={r:t for r,t in terms.items() if t!=0}
    coh=all(t*(-1)**(sum(k)-1)>0 for t in nz.values())
    edge=(-1)**sum(k)*G.H(k,Fr(-Ms-2,2))
    caseB=G.caseB(k)
    st[(caseB,cs,coh,len(nz)>0,edge>0)]+=1
    if not coh or edge<=0: ex.append((k,caseB,{G.nm(r):str(t) for r,t in nz.items()},str(edge)))
print('(caseB, Pieri signs alternate, edge terms coherent, some nonzero term, (-1)^|k| H(-M0/2)>0):')
for kk,v in sorted(st.items()): print(kk,v)
for e in ex[:30]: print(e)
