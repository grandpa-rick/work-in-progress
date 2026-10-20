"""Twisted coherence: for each sign vector s in {+1,-1}^N, is P(s*w) one-signed coefficientwise
(P = Weyl-stripped residual numerator)? Twist s is relative to the pair: we record it as
(sign on w_a, signs on the other free variables in increasing order). Reports, per c, how many incoherent cases are
repaired by SOME twist, and whether a single uniform twist works for all cases (the 'canonical' salvage)."""
import sys,glob,ast,re,itertools
N=int(sys.argv[1])
def parse(P):
    # returns list of (coef, exps) from flint str
    P=P.replace(' ','').replace('-','+-')
    out=[]
    for t in P.split('+'):
        if not t: continue
        coef=1; e=[0]*N
        for f in t.split('*'):
            if f.startswith('-'): coef=-coef; f=f[1:]
            if not f: continue
            if f[0]=='w':
                v,_,p=f.partition('^'); e[int(v[1:])-1]+=int(p) if p else 1
            else: coef*=int(f) if '/' not in f else float(eval(f))
        out.append((coef,e))
    return out
for fn in sorted(glob.glob('coh%d_c*.txt'%N),key=lambda f:int(re.search(r'_c(-?\d+)',f).group(1))):
    c=int(re.search(r'_c(-?\d+)',fn).group(1)); L=open(fn).read()
    if '#' not in L: continue
    inc=0; rep=0; uni=None
    for l in L.splitlines():
        if l.startswith('#'): continue
        key,d=ast.literal_eval(l)
        if d['zero'] or d['p1s']: continue
        inc+=1; a,b=key[0],key[1]; free=[a]+[i for i in range(1,N+1) if i not in (a,b)]
        terms=parse(d['P']); ok=set()
        for s in itertools.product((1,-1),repeat=N-1):
            sg={ (co*eval('*'.join(['1']+[str(s[j]**e[free[j]-1]) for j in range(N-1)])))>0 for co,e in terms}
            if len(sg)==1: ok.add(s)
        if ok: rep+=1
        uni=ok if uni is None else uni&ok
    print('N %d c %3d incoherent %3d  repaired by some twist %3d  uniform twists working for all: %s'%(N,c,inc,rep,sorted(uni) if uni is not None else '-'))
