import sys; sys.path.insert(0,"../prove248"); sys.path.insert(0,'../prove247')
import scert7 as S
S.NV=8
for lim,lam,th in [(0,(1,3),(0,1)),(0,(2,4),(0,1)),(0,(1,3),(1,2)),(0,(1,2,4),(0,0,1)),(0,(1,3,4),(0,1,1)),(0,(1,2),(1,2)),(1,(1,2),(1,0)),(0,(1,3),(2,2)),(1,(1,3),(2,1)),(1,(2,3),(0,-1)),(1,(1,2),(-1,-2)),(1,(2,),(-1,))]:
    out=S.psi_opts(lam,th,7,lim)
    ords=[o[0] for o in out]; ext=min(ords) if lim==0 else max(ords)
    print('lim',lim,lam,th)
    for o in out:
        print('   ','*' if o[0]==ext else ' ',o[0],'child',o[4],o[5])
