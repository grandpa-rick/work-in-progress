"""Post-process cohN_c*.txt: P4 = does the Weyl-stripped residual P take both signs on the positive orthant?
(robust incoherence: no re-orientation of root factors can repair it). Also per-triple tracks across c."""
import sys,glob,ast,random,re
from fractions import Fraction
N=int(sys.argv[1]); R=random.Random(5)
pts=[[Fraction(R.randint(1,60),R.randint(1,60)) for _ in range(N)] for _ in range(300)]
def ev(P,p):
    env={'w%d'%(i+1):p[i] for i in range(N)}
    return eval(P.replace('^','**'),{},env)
res={}
for fn in sorted(glob.glob('coh%d_c*.txt'%N),key=lambda f:int(re.search(r'_c(-?\d+)',f).group(1))):
    c=int(re.search(r'_c(-?\d+)',fn).group(1)); lines=[l for l in open(fn) if not l.startswith('#')]
    if not any(l.startswith('#') for l in open(fn)): print('c',c,'INCOMPLETE',len(lines)); continue
    z=coh=inc=rob=0; rows={}
    for l in lines:
        key,d=ast.literal_eval(l)
        if d['zero']: z+=1; rows[key]='0'; continue
        if d['p1s']: coh+=1; rows[key]='+'; continue
        inc+=1; s=set()
        for p in pts:
            v=ev(d['P'],p); s.add((v>0)-(v<0))
        s.discard(0)
        if len(s)>1: rob+=1; rows[key]='X'
        else: rows[key]='x'
    res[c]=rows
    print('N %d c %3d  coherent %3d  incoherent %3d (robust: P changes sign on w>0: %3d)  zero %3d'%(N,c,coh,inc,rob,z))
if len(sys.argv)>2:
    for key in ast.literal_eval(sys.argv[2]):
        print(key,' '.join('%d:%s'%(c,res[c].get(key,'?')) for c in res))
