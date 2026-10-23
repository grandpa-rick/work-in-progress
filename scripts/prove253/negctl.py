import sys, itertools
sys.path.insert(0,'../wake252/caseB'); sys.path.insert(0,'../wake251/orderpoly')
from brpieri import pieri_dim, boxes
from strip import sp, r, frombox
for n in range(2,int(sys.argv[1])+1):
    agree=tot=zero=0
    for k in sp(n):
        R=r(k)
        if R%2==0 or R<3: continue
        B=boxes(k); th={x for x in B if x[0]+x[1]>=R-2}; rs=frombox(B-th)
        conflict=any(a==b and a+b==R-1 and ((a+1,b) in th or (a,b+1) in th) for (a,b) in th)
        # partners of outer box (a,a): inner boxes in same row/col: (a,a-1) invalid, (a-1,a) above
        conflict=any(a==b and a+b==R-1 and (a-1,a) in th for (a,b) in th)
        c=pieri_dim(k,n).get(rs,0); tot+=1; zero+= c==0
        agree += (c==0)==conflict
    print(f'm={n}: r odd, 2-strip: predicted-zero matches BR: {agree}/{tot}, zeros {zero}')
