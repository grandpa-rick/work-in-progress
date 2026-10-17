"""Assemble Ask-1 table: actual full-q poles (polesnum pn_*.txt) vs Clio's predicted pairs, stratum status per label,
q-valuation of the residue from day243/q488c/qinf_N5_partial.txt (symbolic q, first-order residue)."""
import glob, ast, re, itertools, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../day243/q488c'))
from engine import subsets, geq
D=os.path.dirname(os.path.abspath(__file__)); N=5; S=subsets(N)
def rm(M,a,b): return tuple(x for x in M if x not in (a,b))
def pred(l,m,nu):
    return {(a,b) for k in S if geq(nu,k) for a,b in itertools.combinations(k,2) if geq(rm(k,a,b),l) and geq(rm(k,a,b),m)}
qinf={}
for line in open(os.path.join(D,'../../day243/q488c/qinf_N5_partial.txt')):
    i=line.index('[('); head=line[:i].replace(') (',')|(').split('|')
    t=tuple(ast.literal_eval(x.strip()) for x in head)
    t=tuple(tuple(x) if isinstance(x,tuple) else (x,) for x in t)
    qinf[t]={p:v for p,v,_ in ast.literal_eval(line[i:])}
P={}
for f in glob.glob(os.path.join(D,'pn_*.txt')):
    for line in open(f):
        mm=re.match(r'(\d+) (.*) \((\d), (\d)\) pole order .*: \[(\w+), (\w+)\]',line)
        if not mm: continue
        idx=int(mm.group(1)); a,b=int(mm.group(3)),int(mm.group(4)); o=(mm.group(5),mm.group(6))
        P.setdefault(idx,{})[(a,b)]=o
T=ast.literal_eval(open(os.path.join(D,'../../day243/q488c/n5fails_partial.txt')).read())
mx=lambda X: max(X) if X else 0
agree=True; npole=0; inS={'lam':0,'mu':0,'both':0}
for i,(l,m,nu) in enumerate(T):
    pr=pred(l,m,nu); got=P.get(i,{})
    if len(got)<10: print(i,'INCOMPLETE',len(got)); agree=False; continue
    poles={p:o for p,o in got.items() if o!=('0','0')}
    if set(poles)!=pr: agree=False
    rows=[]
    for (a,b),o in sorted(poles.items()):
        npole+=1
        il=mx(l)<=b; im=mx(m)<=b
        inS['lam']+=il; inS['mu']+=im; inS['both']+=(il and im)
        tq=tuple(tuple(x) for x in (l,m,nu))
        rows.append('H%d%d ord %s/%s  I=lam:%s I=mu:%s  res qval0 %s'%(a,b,o[0],o[1],'in' if il else 'OUT','in' if im else 'OUT',qinf.get(tq,{}).get((a,b),'?')))
    print(l,m,nu,'maxlam',mx(l),'maxmu',mx(m),'poles==predicted:',set(poles)==pr)
    for r in rows: print('   ',r)
print('total pole (triple,pair) instances',npole,'| in-stratum for I=lam',inS['lam'],'I=mu',inS['mu'],'both',inS['both'])
print('actual pole set == Clio Conj divisor predicted set on all 12:',agree)
