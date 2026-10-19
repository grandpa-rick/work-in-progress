import scert
from scert import *
orig=scert.certify
def certify(S,n,a,b,cr,depth=0):
    S2=scert.clean(S,n,a,b)
    r=orig(S,n,a,b,cr,depth)
    print('  '*depth+f'n={n} cr={cr} #S={len(S2)} -> {r}', [ (v,k[1],k[2],k[3],k[4],k[5],k[6],k[7]) for k,v in list(S2.items())[:4]])
    return r
scert.certify=certify
start('C',5,(),(1,2),(4,5),None,None)
