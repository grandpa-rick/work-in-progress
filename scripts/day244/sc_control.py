from sc import *
B=basis('K',(1,2,0))   # lam=3/4: alpha=0, printed (typo) 23|e13
R=sc(B); print('Koncki lam=3/4 PRINTED:',{f for c in R.values() for e in c.values() for f in poles(e)})
B['23','13']=-(y+1)*(s*t-1)*(s**2*t*y+s*y+s+t)/t   # LRR value, t^{1-ceil2}=t^{-1} times LRR bracket
R=sc(B); print('Koncki lam=3/4 LRR-fixed:',{f for c in R.values() for e in c.values() for f in poles(e)})
