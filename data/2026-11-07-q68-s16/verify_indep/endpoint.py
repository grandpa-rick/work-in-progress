import sys,json
from vi import *
def why(a,b):
    la,lb=length(a),length(b); lo,hi=(a,b) if la<lb else (b,a)
    if la==lb: return 'equal length'
    if not bruhat_le(lo,hi): return 'incomparable (gap %d)'%abs(la-lb)
    if abs(la-lb)%2==0: return 'even gap %d'%abs(la-lb)
    if abs(la-lb)==1: return 'gap1 comparable (mu=1)'
    return 'F1 violation gap %d'%abs(la-lb)
for f in sys.argv[1:]:
    for r in json.load(open(f)):
        if not r['parity']: continue
        x=inv_rsk(r['Px'],r['Q']); z=inv_rsk(r['Pz'],r['Q'])
        for s in r['kills'][:4]:
            y=lmul(s,x); m,path=star_search(y,z)
            a,b=y,z
            for j in path: a,b=lstar(j,a),lstar(j,b)
            print(f, 'Q=',r['Q'][:2],'...','s=',s,'l(sx),l(z)=',length(y),length(z),'path',path,'->',why(a,b))
